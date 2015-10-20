import pysam
import vcf
import numpy as np
import csv
import ipdb
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

    split_reads = np.where([process.is_supporting_split_read_lenient(x,pos,params.tr*2) for x in loc_reads])[0]
    split_all = loc_reads[split_reads]
    
    if len(split_all)>0:
        pos, total = sum(split_all['align_start'] < sc_len), len(split_all)
        if pos/float(total) > 0.1 and pos/float(total) < 0.9:
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

def get_dir(loc_reads,pos):
    split_reads = np.where([process.is_supporting_split_read_lenient(x,pos) for x in loc_reads])[0]
    split_all = loc_reads[split_reads]
    if len(split_all)>0:
        dir_split = get_dir_split(split_all,params.tr)
        return dir_split
    else:
        return '?'

def get_consensus_align(loc_reads,pos):
    split_reads = np.where([process.is_supporting_split_read_lenient(x,pos,params.tr*2) for x in loc_reads])[0]
    split_all = loc_reads[split_reads]
   
    if len(split_all)!=0:
        ref_starts = [x for x in split_all['ref_start']]
        ref_ends   = [x for x in split_all['ref_end']  ]
        
        consensus_align_right = max(set(ref_ends), key=ref_ends.count)
        consensus_align_left  = max(set(ref_starts), key=ref_starts.count)

        return consensus_align_right, consensus_align_left+1
    else:
        return 0,0

def get_bp_dir(sv,loc_reads,pos,sc_len,bp_num):

    bp_dir = 'bp%d_dir' % bp_num
    sv_class = str(sv['classification'])
    
    ca_right, ca_left = get_consensus_align(loc_reads,pos)
    ca_right = ca_right if (ca_right-params.tr*2 < pos and ca_right+params.tr*2 > pos) else 0
    ca_left  = ca_left  if (ca_left-params.tr*2 < pos  and ca_left+params.tr*2 > pos)  else 0

    mixed_right = has_mixed_evidence(loc_reads,ca_right,sc_len) if ca_right != 0 else False
    mixed_left = has_mixed_evidence(loc_reads,ca_right,sc_len) if ca_right != 0 else False
    
    if has_mixed_evidence(loc_reads,pos,sc_len) or mixed_right or mixed_left:
        sv[bp_dir] = '?'
        sv['classification'] = 'MIXED' if sv_class=='' else sv_class+';MIXED'
    else:
        sv[bp_dir] = get_dir(loc_reads,pos)
        if sv[bp_dir] == '?':
            sv['classification'] = 'UNKNOWN_DIR' if sv_class=='' else sv_class+';UNKNOWN_DIR'

    return sv, ca_right, ca_left

def get_dir_info(row,bam,max_dep):

    sv = np.array(row,copy=True)
    bp_dtype = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S1')]
    
    sv_id, bp1_chr, bp1_pos, bp1_dir, bp2_chr, bp2_pos, bp2_dir, sv_class = [h[0] for h in params.sv_dtype]
    bp1 = np.array((sv[bp1_chr],sv[bp1_pos]-(params.tr*2),sv[bp1_pos]+(params.tr*2),sv[bp1_dir]),dtype=bp_dtype)
    bp2 = np.array((sv[bp2_chr],sv[bp2_pos]-(params.tr*2),sv[bp2_pos]+(params.tr*2),sv[bp2_dir]),dtype=bp_dtype)
    
    bamf = pysam.AlignmentFile(bam, "rb")
    loc1_reads, err_code1 = process.get_loc_reads(bp1,bamf,max_dep)    
    loc2_reads, err_code2 = process.get_loc_reads(bp2,bamf,max_dep)    
    bamf.close()

    sv_class = str(sv['classification'])
    if err_code1 == 1 or err_code2 == 1:
        sv['classification'] = 'HIDEP' if sv_class=='' else sv_class+';HIDEP' 
        return sv, (0,0,0,0)
    elif err_code1 == 2 or err_code2 == 2:
        sv['classification'] = 'READ_FETCH_FAILED' if sv_class=='' else sv_class+';READ_FETCH_FAILED'
        return sv, (0,0,0,0)

    sv, bp1_ca_right, bp1_ca_left = get_bp_dir(sv,loc1_reads,sv['bp1_pos'],params.tr,1)
    sv, bp2_ca_right, bp2_ca_left = get_bp_dir(sv,loc2_reads,sv['bp2_pos'],params.tr,2)
    consens_aligns = (bp1_ca_right,bp1_ca_left,bp2_ca_right,bp2_ca_left)

    return sv, consens_aligns

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

def set_dir_class(sv,bp1_dir,bp2_dir,sv_class,new_pos1,new_pos2):
    sv['bp1_dir'] = bp1_dir
    sv['bp2_dir'] = bp2_dir
    sv['classification'] = sv_class
    if new_pos1!=0:
        sv['bp1_pos'] = new_pos1
    if new_pos2!=0:
        sv['bp2_pos'] = new_pos2
    return sv

def is_same_sv(sv1,sv2,threshold=params.tr):
    sv1_chr1, sv1_bp1, sv1_dir1, sv1_chr2, sv1_bp2, sv1_dir2 = sv1
    sv2_chr1, sv2_bp1, sv2_dir1, sv2_chr2, sv2_bp2, sv2_dir2 = sv2

    if sv1_chr1==sv2_chr1 and sv1_chr2==sv2_chr2 and sv1_dir1 == sv2_dir1 and sv2_dir2 == sv2_dir2:
        if abs(sv1_bp1-sv2_bp1)<threshold and abs(sv1_bp2-sv2_bp2)<threshold:
            return True
    if sv1_chr2==sv2_chr1 and sv1_chr1==sv2_chr2 and sv1_dir2 == sv2_dir1 and sv2_dir1 == sv2_dir2:
        if abs(sv1_bp2-sv2_bp1)<threshold and abs(sv1_bp1-sv2_bp2)<threshold:
            return True
    return False

def remove_duplicates(svs):
    to_delete = []    
    for idx1,sv1 in enumerate(svs):
        if idx1+1 == len(svs): break
        for idx2,sv2 in enumerate(svs[idx1+1:]):
            if idx1!=idx2:
                sv1_tmp = (sv1['bp1_chr'],sv1['bp1_pos'],sv1['bp1_dir'],\
                           sv1['bp2_chr'],sv1['bp2_pos'],sv1['bp2_dir'])
                sv2_tmp = (sv2['bp1_chr'],sv2['bp1_pos'],sv2['bp1_dir'],\
                           sv2['bp2_chr'],sv2['bp2_pos'],sv1['bp2_dir'])
                if is_same_sv(sv1_tmp,sv2_tmp):
                    to_delete.append(idx1)
    svs = np.delete(svs,np.array(to_delete))
    return svs

def split_mixed_svs(svs,ca):
    for idx,row in enumerate(svs):
        sv_class = np.array(row['classification'].split(';'))
        if 'MIXED' in sv_class: 
            if 'UNKNOWN_DIR' in sv_class:
                svs[idx]['classification'] = 'UNKNOWN_DIR'
                continue #can't fix this if one direction is unknown
            if len(sv_class)>1 and sv_class[0]=='MIXED' and sv_class[1]=='MIXED':
                #likely an inversion
                #set dirs to + and create new sv with - dirs 
                new_sv = svs[idx].copy()
                new_pos1, new_pos2 = ca[idx]['bp1_ca_right'],ca[idx]['bp2_ca_right']
                new_sv = set_dir_class(new_sv,'+','+','',new_pos1, new_pos2)                
                new_sv['ID'] = svs[len(svs)-1]['ID']+1
                svs = np.append(svs,new_sv)

                new_pos1, new_pos2 = ca[idx]['bp1_ca_left'],ca[idx]['bp2_ca_left']
                svs[idx] = set_dir_class(svs[idx],'-','-','',new_pos1,new_pos2)                

            elif row['bp1_dir']!='?':
                #set dir to -, create new sv with dir set to +
                new_sv = svs[idx].copy()
                new_pos2 = ca[idx]['bp2_ca_right']
                new_sv = set_dir_class(new_sv,row['bp1_dir'],'+','',0,new_pos2)
                svs = np.append(svs,new_sv)

                new_pos2 = ca[idx]['bp2_ca_left']
                svs[idx] = set_dir_class(svs[idx],row['bp1_dir'],'-','',0,new_pos2)

            elif row['bp2_dir']!='?':
                #set dir to -, create new sv with dir set to +
                new_sv = svs[idx].copy()
                new_pos1 = ca[idx]['bp1_ca_right']
                new_sv = set_dir_class(new_sv,'+',row['bp2_dir'],'',new_pos1,0)
                svs = np.append(svs,new_sv)
                
                svs[idx] = set_dir_class(svs[idx],'-',row['bp2_dir'],'',new_pos1,0)

            else:
                svs[idx]['classification'] = 'UNKNOWN_DIR'
    
    svs = remove_duplicates(svs)
    return svs

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
   
    consens_dtype = [('bp1_ca_right',int),('bp1_ca_left',int),('bp2_ca_right',int),('bp2_ca_left',int)]
    ca = np.zeros(len(svs),dtype=consens_dtype)
        
    if not use_dir:
        for idx,sv in enumerate(svs):
            svs[idx], ca[idx]  = get_dir_info(sv,bam,max_dep)            
        svs = split_mixed_svs(svs,ca)

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
        header = [field for field,dtype in params.sv_dtype]
        writer.writerow(header)
        for sv_out in svs:
            writer.writerow(sv_out)
