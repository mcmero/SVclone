import pysam
import vcf
import numpy as np
import ipdb
import csv
from . import load_data
from . import process
from . import parameters as params
from . import svDetectFuncs as svd

def classify_event(sv,sv_id,svd_prevResult,prevSV):    
    
    svd_result = svd.detect(prevSV,svd_prevResult,sv)
    svd_prevResult,prevSV = svd_result,sv

    classification = svd.getResultType(svd_result)
    sv_id = sv_id if svd_result[0]==svd.SVtypes.interspersedDuplication else sv_id+1
    
    return classification, sv_id, svd_prevResult, prevSV

def get_dir_split(split,sc_len):
    align_mean =  np.mean(split['align_start'])
    assign_dir = '+' if align_mean < sc_len else '-'    
    return assign_dir

#TODO: optimise calling direction
def get_dir(loc_reads,pos):
    split_reads = np.where([process.is_supporting_split_read_lenient(x,pos) for x in loc_reads])[0]
    split_all = loc_reads[split_reads]
    if len(split_all)>0:
        dir_split = get_dir_split(split_all,params.tr)
        return dir_split
    else:
        return '?'

def get_dir_info(row,bam,max_dep):

    sv = np.array(row,copy=True)
    bp_dtype = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S1')]
    
    sv_id, bp1_chr, bp1_pos, bp1_dir, bp2_chr, bp2_pos, bp2_dir, sv_class = [h[0] for h in params.sv_dtype]
    bp1 = np.array((sv[bp1_chr],sv[bp1_pos]-params.tr,sv[bp1_pos]+params.tr,sv[bp1_dir]),dtype=bp_dtype)
    bp2 = np.array((sv[bp2_chr],sv[bp2_pos]-params.tr,sv[bp2_pos]+params.tr,sv[bp2_dir]),dtype=bp_dtype)
    
    bamf = pysam.AlignmentFile(bam, "rb")
    loc1_reads, err_code1 = process.get_loc_reads(bp1,bamf,max_dep)    
    loc2_reads, err_code2 = process.get_loc_reads(bp2,bamf,max_dep)    
    bamf.close()

    if err_code1 == 1 or err_code2 == 1:
        sv['classification'] = 'HIDEP'
        return sv
    elif err_code1 == 2 or err_code2 == 2:
        sv['classification'] = 'READ_FETCH_FAILED'

    sv['bp1_dir'] = get_dir(loc1_reads,sv[bp1_pos])
    sv['bp2_dir'] = get_dir(loc2_reads,sv[bp2_pos])

    return sv
    
def preproc_svs(args):
    
    svin         = args.svin
    bam          = args.bam
    out          = args.out
    max_dep      = args.max_dep
    simple       = args.simple_svs
    socrates     = args.socrates
    use_dir      = args.use_dir
    min_mapq     = args.min_mapq
    class_field  = args.class_field
    filt_repeats = args.filt_repeats
    
    filt_repeats = filt_repeats.split(',') if filt_repeats != '' else filt_repeats
    filt_repeats = [rep for rep in filt_repeats if rep!='']    
    
    outf = '%s_svin.txt'%out
    
    svs = np.empty(0)
    if simple:
        svs = load_data.load_input_simple(svin,use_dir,class_field)
    elif socrates:
        svs = load_data.load_input_socrates(svin,use_dir,min_mapq,filt_repeats)
    else:
        svs = load_data.load_input_vcf(svin,class_field)

    sv_id = 0
    svd_prevResult, prevSV = None, None
    sv_input = np.empty(0,dtype=params.sv_dtype)
    # add directionality information
    for row in svs:
        sv = get_dir_info(row,bam,max_dep)
        if class_field!='' and sv['classification']!='':
            sv_id += 1
        if sv['classification']=='' and (sv['bp1_dir']=='?' or sv['bp2_dir']=='?'):
            sv_id += 1
            sv['classification'] = 'UNKNOWN_DIR'
        elif sv['classification']=='':
            sv['classification'], sv_id, svd_prevResult, prevSV = classify_event(sv,sv_id,svd_prevResult,prevSV)
        else:
            sv_id += 1
        sv['ID'] = sv_id
        sv_input = np.append(sv_input,sv)
        
    # reclassify once all have been processed
    for sv in sv_input:
        trx_label    = svd.getResultType([svd.SVtypes.translocation])
        intins_label = svd.getResultType([svd.SVtypes.interspersedDuplication])
        if sv['classification']==intins_label:
            sv_info[idx-1]['classification'] = intins_label
            translocs = svd.detectTransloc(idx,sv_info)
            if len(translocs)>0:
                for i in translocs: sv_input[i]['classification'] = trx_label

    with open('%s_svin.txt'%out,'w') as outf:
        writer = csv.writer(outf,delimiter='\t',quoting=csv.QUOTE_NONE,quotechar='')
        header = ['ID','bp1_chr','bp1_pos','bp1_dir','bp2_chr','bp2_pos', 'bp2_dir', 'classification']       
        writer.writerow(header)
        for sv_out in sv_input:
            writer.writerow(sv_out)
    
