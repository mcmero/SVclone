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

def has_mixed_evidence(loc_reads,pos,sc_len):
    split_reads = np.where([process.is_supporting_split_read_lenient(x,pos) for x in loc_reads])[0]
    split_all = loc_reads[split_reads]
    
    if len(split_all)>0:
        pos, total = sum(split_all['align_start'] < sc_len), len(split_all)
        if pos/float(total) > 0.2 and pos/float(total) < 0.8:
            return True

    return False

def get_dir_span(span):
    is_rev = np.sum(span['is_reverse'])
    assign_dir = '+' if is_rev <= len(span)/2 else '-'
    return assign_dir

def get_dir_split(split,sc_len):
    #align_mean =  np.mean(split['align_start'])    
    #assign_dir = '+' if align_mean < sc_len else '-'
    pos, total = sum(split['align_start'] < sc_len), len(split)
    assign_dir = '+' if pos/float(total) > 0.5 else '-'
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

def get_bp_dir(sv,loc_reads,pos,sc_len,bp_num):

    bp_dir = 'bp%d_dir' % bp_num
    sv_class = str(sv['classification'])

    if has_mixed_evidence(loc_reads,pos,sc_len):
        sv[bp_dir] = '?'
        sv['classification'] = 'MIXED' if sv_class=='' else sv_class+';MIXED'
    else:
        sv[bp_dir] = get_dir(loc_reads,pos)
        if sv[bp_dir] == '?':
            sv['classification'] = 'UNKNOWN_DIR' if sv_class=='' else sv_class+';UNKNOWN_DIR'

    return sv

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

    sv_class = str(sv['classification'])
    if err_code1 == 1 or err_code2 == 1:
        sv['classification'] = 'HIDEP' if sv_class=='' else sv_class+';HIDEP' 
        return sv
    elif err_code1 == 2 or err_code2 == 2:
        sv['classification'] = 'READ_FETCH_FAILED' if sv_class=='' else sv_class+';READ_FETCH_FAILED'
        return sv

    sv = get_bp_dir(sv,loc1_reads,sv['bp1_pos'],params.tr,1)
    sv = get_bp_dir(sv,loc2_reads,sv['bp2_pos'],params.tr,2)

    return sv
    
#def mixed_svs_remain(svs):
#    sv_classes = map(lambda x: x.split(';'),svs['classification'])
#    sv_classes [svc for svc in sv_classes]
#    return 'MIXED' in sv_classes:
#
#def does_break_match(chr1,pos1,chr2,pos2):
#    return chr1==chr2 and (pos1 + params.tr) > pos2 and (pos1 - params.tr) < pos2:
#
#def get_matching_svs(bp_chr,bp_pos,svs):
#    sv_matches = np.empty(0,dtype=params.sv_dtype)
#    which_matches = []
#    for idx,row in enumerate(svs):
#        matches = []
#        if not np.all(row==sv) and does_break_match(bp_chr,bp_pos,row['bp1_chr'],row['bp1_pos']):
#            sv_matches = np.append(sv_matches,row)
#            matches.append([idx,1])
#        if not np.all(row==sv) and does_break_match(bp_chr,bp_pos,row['bp2_chr'],row['bp2_pos']):
#            sv_matches = np.append(sv_matches,row)
#            matches.append([idx,2])
#        which_matches.append(matches)
#    return sv_matches,which_matches

def set_dir_class(sv,bp1_dir,bp2_dir,sv_class):
    sv['bp1_dir'] = bp1_dir
    sv['bp2_dir'] = bp2_dir
    sv['classification'] = sv_class
    return sv

def split_mixed_svs(svs):
    for idx,row in enumerate(svs):
        sv_class = np.array(row['classification'].split(';'))
        if 'MIXED' in sv_class:
            if 'UNKNOWN_DIR' in sv_class:
                svs[idx]['classification'] = 'UNKNOWN_DIR'
                continue #can't fix this if one direction is unknown
            if len(sv_class)>1 and sv_class[0]=='MIXED' and sv_class[1]=='MIXED':
                #likely an inversion
                #set dirs to + and create new sv with - dirs 
                svs[idx] = set_dir_class(svs[idx],'-','-','')
                new_sv = np.array(svs[idx],copy=True)
                new_sv = set_dir_class(new_sv,'+','+','')
                svs = np.append(svs,new_sv)
            elif row['bp1_dir']!='?':
                #set dir to -, create new sv with dir set to +
                svs[idx]['bp2_dir'] = '-'
                new_sv = np.array(svs[idx],copy=True)
                new_sv = set_dir_class(new_sv,row['bp1_dir'],'+','')
                svs = np.append(svs,new_sv)
            elif row['bp2_dir']!='?':
                #set dir to -, create new sv with dir set to +
                svs[idx]['bp1_dir'] = '-'
                new_sv = np.array(svs[idx],copy=True)
                new_sv = set_dir_class(new_sv,'+',row['bp2_dir'],'')
                svs = np.append(svs,new_sv)
    return np.unique(svs)

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
    
    if not use_dir:
        for idx,sv in enumerate(svs):
            svs[idx] = get_dir_info(sv,bam,max_dep)
        svs = split_mixed_svs(svs)

    sv_id = 0
    svd_prevResult, prevSV = None, None
    for idx,sv in enumerate(svs):
        if sv['classification']=='':
            sv_class, sv_id, svd_prevResult, prevSV = classify_event(sv,sv_id,svd_prevResult,prevSV)
            svs[idx]['classification'] = sv_class
        else:
            sv_id += 1
        svs[idx]['ID'] = sv_id
    
    # reclassify once all have been processed
    for idx,sv in enumerate(svs):
        trx_label    = svd.getResultType([svd.SVtypes.translocation])
        intins_label = svd.getResultType([svd.SVtypes.interspersedDuplication])
        if sv['classification']==intins_label:
            svs[idx-1]['classification'] = intins_label
            translocs = svd.detectTransloc(idx,svs)
            if len(translocs)>0:
                for i in translocs: svs[i]['classification'] = trx_label

    with open('%s_svin.txt'%out,'w') as outf:
        writer = csv.writer(outf,delimiter='\t',quoting=csv.QUOTE_NONE,quotechar='')
        header = ['ID','bp1_chr','bp1_pos','bp1_dir','bp2_chr','bp2_pos', 'bp2_dir', 'classification']       
        writer.writerow(header)
        for sv_out in svs:
            writer.writerow(sv_out)
