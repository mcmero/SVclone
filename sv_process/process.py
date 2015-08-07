'''
Using characterised SVs, count normal and supporting reads at SV locations
'''

import os
import numpy as np
import pysam
import csv
from collections import OrderedDict
from operator import methodcaller
import vcf
from . import parameters as params
from . import bamtools
from . import svDetectFuncs as svd

#TODO: remove after testing
import ipdb

def read_to_array(x,bamf):
    chrom = bamf.getrname(x.reference_id)
    try:
        read = np.array((x.query_name,chrom,x.reference_start,x.reference_end,x.query_alignment_start,
                         x.query_alignment_end,x.query_length,x.tlen,np.bool(x.is_reverse)),dtype=params.read_dtype)
        return read
    except TypeError:
        print 'Warning: record %s contains invalid attributes' % x.query_name
        return np.empty(len(params.read_dtype),dtype=params.read_dtype)

def is_soft_clipped(r):
    return r['align_start'] != 0 or (r['len'] + r['ref_start'] != r['ref_end'])

def is_minor_softclip(r):
    return (r['align_start'] < (params.tr)) and ((r['align_end'] + r['ref_start'] - r['ref_end']) < (params.tr))

def is_normal_across_break(r,pos,max_ins,sc_len):
    # must overhang break by at least soft-clip threshold
    return  (not is_soft_clipped(r)) and \
            (abs(r['ins_len']) < max_ins) and \
            (r['ref_start'] <= (pos - sc_len)) and \
            (r['ref_end'] >= (pos + sc_len)) 
          
def is_normal_spanning(r,m,pos,max_ins):
    if not (is_soft_clipped(r) or is_soft_clipped(m)):
        if (not r['is_reverse'] and m['is_reverse']) or (r['is_reverse'] and not m['is_reverse']):
            return (abs(r['ins_len']) < max_ins) and \
                   (r['ref_end'] < (pos + params.tr)) and \
                   (m['ref_start'] > (pos - params.tr))
    return False

def is_supporting_split_read(r,pos,max_ins,sc_len):
    '''
    Return whether read is a supporting split read.
    Doesn't yet check whether the soft-clip aligns
    to the other side.
    '''
    if r['align_start'] < (params.tr): #a "soft" threshold if it is soft-clipped at the other end        
        return r['ref_end'] > (pos - params.tr) and r['ref_end'] < (pos + params.tr) and \
            (r['len'] - r['align_end'] >= sc_len) and abs(r['ins_len']) < max_ins
    else:
        return r['ref_start'] > (pos - params.tr) and r['ref_start'] < (pos + params.tr) and \
            (r['align_start'] >= sc_len) and abs(r['ins_len']) < max_ins

def is_supporting_split_read_lenient(r,pos):
    '''
    Same as is_supporting_split_read wihout insert and soft-clip threshold checks
    '''
    if r['align_start'] < (params.tr): #a "soft" threshold if it is soft-clipped at the other end        
        return (r['len'] - r['align_end'] >= params.tr) and r['ref_end'] > (pos - params.tr) and r['ref_end'] < (pos + params.tr)
    else:
        return (r['align_start'] >= params.tr) and r['ref_start'] > (pos - params.tr) and r['ref_start'] < (pos + params.tr)

def get_sc_bases(r,pos):
    '''
    Return the number of soft-clipped bases
    '''
    if r['align_start'] < (params.tr):
        return r['len'] - r['align_end']
    else:
        return r['align_start']

def get_bp_dist(x,bp_pos):
    if x['is_reverse']: 
        return (x['ref_end'] - bp_pos)
    else: 
        return (bp_pos - x['ref_start'])

def is_supporting_spanning_pair(r,m,bp1,bp2,inserts,max_ins):
    pos1 = (bp1['start'] + bp1['end']) / 2
    pos2 = (bp2['start'] + bp2['end']) / 2
    
    if is_soft_clipped(r) or is_soft_clipped(m):
        return False
    
    #ensure this isn't just a regular old spanning pair    
    if r['chrom']==m['chrom']:
        if r['ref_start']<m['ref_start']:
            if m['ref_start']-r['ref_end'] < max_ins: return False
        else:
            if r['ref_start']-m['ref_end'] < max_ins: return False

    ins_dist1 = get_bp_dist(r,pos1)
    ins_dist2 = get_bp_dist(m,pos2)

    if ins_dist1<0 or ins_dist2<0:
        return False
    else:
        return (ins_dist1+ins_dist2) < max_ins

def get_loc_reads(bp,bamf,max_dp):
    loc = '%s:%d:%d' % (bp['chrom'], max(0,bp['start']), bp['end'])
    loc_reads = np.empty([0,len(params.read_dtype)],dtype=params.read_dtype)    
    try:
        iter_loc = bamf.fetch(region=loc,until_eof=True)
        for x in iter_loc:
            read = read_to_array(x,bamf) 
            loc_reads = np.append(loc_reads,read)
            if len(loc_reads) > max_dp:
                print('Read depth too high at %s' % loc)
                return np.empty(0)

        loc_reads = np.sort(loc_reads,axis=0,order=['query_name','ref_start'])
        loc_reads = np.unique(loc_reads) #remove duplicates
        return loc_reads
    except ValueError:
        print('Fetching reads failed for loc: %s' % loc)
        return np.empty(0)

def reads_to_sam(reads,bam,bp1,bp2,name):
    '''
    For testing read assignemnts.
    Takes reads from array, matches them to bam 
    file reads by query name and outputs them to Sam
    '''
    if len(reads)==0:
        return None

    bamf = pysam.AlignmentFile(bam, "rb")
    loc1 = '%s:%d:%d' % (bp1['chrom'], bp1['start'], bp1['end'])
    loc2 = '%s:%d:%d' % (bp2['chrom'], bp2['start'], bp2['end'])
    iter_loc1 = bamf.fetch(region=loc1,until_eof=True)
    iter_loc2 = bamf.fetch(region=loc2,until_eof=True)
    
    loc1 = '%s-%d' % (bp1['chrom'], (bp1['start']+bp1['end'])/2)
    loc2 = '%s-%d' % (bp2['chrom'], (bp1['start']+bp1['end'])/2)
    sam_name = '%s_%s-%s' % (name,loc1,loc2)
    bam_out = pysam.AlignmentFile('%s.sam'%sam_name, "w", header=bamf.header)
    
    for x in iter_loc1:
        [bam_out.write(x) for r in reads if r['query_name']==x.query_name]
    for x in iter_loc2:
        [bam_out.write(x) for r in reads if r['query_name']==x.query_name]
    
    bamf.close()
    bam_out.close()

def windowed_norm_read_count(loc_reads,inserts,max_ins):
    '''
    Counts normal non-soft-clipped reads within window range
    '''
    cnorm = 0
    for idx,r in enumerate(loc_reads):
        if idx+1 >= len(loc_reads):
            break    
        r1 = np.array(loc_reads[idx],copy=True)
        r2 = np.array(loc_reads[idx+1],copy=True)
        if r1['query_name']!=r2['query_name'] or r1['chrom']!=r2['chrom']:
            continue
        ins_dist = r2['ref_end']-r1['ref_start']
        facing = not r1['is_reverse'] and r2['is_reverse']
        if not is_soft_clipped(r1) and not is_soft_clipped(r2) and facing and ins_dist > 0 and ins_dist < max_ins:
            cnorm = cnorm + 2
    return cnorm

def get_loc_counts(loc_reads,pos,rc,reproc,split,norm,max_ins,sc_len,bp_num=1):
    for idx,x in enumerate(loc_reads):
        if idx+1 >= len(loc_reads):            
            break        
        r1 = loc_reads[idx]
        r2 = loc_reads[idx+1] if (idx+2)<=len(loc_reads) else None
        
        if is_normal_across_break(x,pos,max_ins,sc_len):
            norm = np.append(norm,r1)            
            split_norm = 'bp%d_split_norm'%bp_num
            rc[split_norm] = rc[split_norm]+1 
        elif is_supporting_split_read(x,pos,max_ins,sc_len):
            split = np.append(split,x)            
            split_supp = 'bp%d_split'%bp_num
            split_cnt = 'bp%d_sc_bases'%bp_num
            rc[split_supp] = rc[split_supp]+1 
            rc[split_cnt] = rc[split_cnt]+get_sc_bases(x,pos)
        elif r2!=None and r1['query_name']==r2['query_name'] and is_normal_spanning(r1,r2,pos,max_ins):
            norm = np.append(norm,r1)            
            norm = np.append(norm,r2)            
            span_norm = 'bp%d_span_norm'%bp_num
            rc[span_norm] = rc[span_norm]+1 
        else:
            reproc = np.append(reproc,x) #may be spanning support or anomalous
    return rc, reproc, split, norm

def get_dir_split(split,sc_len):
    align_mean =  np.mean(split['align_start'])
    assign_dir = '+' if align_mean < sc_len else '-'
    return assign_dir

def get_dir_span(span):
    is_rev = np.sum(span['is_reverse'])
    assign_dir = '+' if is_rev <= len(span)/2 else '-'
    return assign_dir

def get_dir(split,span,loc_reads,pos,sc_len):    
    # split read direction tends to be more reliable
    if len(split)>0 and len(split)>=len(span):
        dir_split = get_dir_split(split,sc_len)
        return dir_split
        
    elif len(span)>0:
        dir_span = get_dir_span(span)
        return dir_span 
        
    else:
        split_reads = np.where([is_supporting_split_read_lenient(x,pos) for x in loc_reads])[0]
        split_all = loc_reads[split_reads]
        if len(split_all)>0:
            dir_split = get_dir_split(split_all,sc_len)
            return dir_split
        else:
            return '?'

def get_sv_read_counts(bp1,bp2,bam,inserts,max_dp,max_ins,sc_len,use_dir):
    bamf = pysam.AlignmentFile(bam, "rb")
    pos1 = (bp1['start'] + bp1['end']) / 2
    pos2 = (bp2['start'] + bp2['end']) / 2
    loc1_reads = get_loc_reads(bp1,bamf,max_dp)    
    loc2_reads = get_loc_reads(bp2,bamf,max_dp)    
    bamf.close()
    if len(loc1_reads)==0 or len(loc2_reads)==0:
        return np.empty(0,dtype=params.sv_out_dtype)
    
    rc = np.zeros(1,dtype=params.sv_out_dtype)
    #rc['sv_id']='%s:%d-%s:%d'%(bp1['chrom'],pos1,bp2['chrom'],pos2)
    reproc = np.empty([0,len(params.read_dtype)],dtype=params.read_dtype)
    
    split_bp1 = np.empty([0,len(params.read_dtype)],dtype=params.read_dtype)
    split_bp2 = np.empty([0,len(params.read_dtype)],dtype=params.read_dtype)
    norm = np.empty([0,len(params.read_dtype)],dtype=params.read_dtype)
    
    rc, reproc, split_bp1, norm = get_loc_counts(loc1_reads,pos1,rc,reproc,split_bp1,norm,max_ins,sc_len)
    rc, reproc, split_bp2, norm = get_loc_counts(loc2_reads,pos2,rc,reproc,split_bp2,norm,max_ins,sc_len,2)
    rc['bp1_win_norm'] = windowed_norm_read_count(loc1_reads,inserts,max_ins)
    rc['bp2_win_norm'] = windowed_norm_read_count(loc2_reads,inserts,max_ins)
    
    reproc = np.sort(reproc,axis=0,order=['query_name','ref_start'])
    reproc = np.unique(reproc) #remove dups
    span_bp1 = np.empty([0,len(params.read_dtype)],dtype=params.read_dtype)
    span_bp2 = np.empty([0,len(params.read_dtype)],dtype=params.read_dtype)
    
    for idx,x in enumerate(reproc):
        if idx+1 >= len(reproc):
            break        
        if reproc[idx+1]['query_name'] != reproc[idx]['query_name']:
            #not paired
            continue
        mate = np.array(reproc[idx+1],copy=True)
        r1 = np.array(x,copy=True)
        r2 = np.array(mate,copy=True)

        #if read corresponds to bp2 and mate to bp1
        if (bp1['chrom']!=bp2['chrom'] and x['chrom']==bp2['chrom']) or \
            (pos1 > pos2 and bp1['chrom']==bp2['chrom']):
            r1 = mate
            r2 = np.array(x,copy=True)
        if is_supporting_spanning_pair(r1,r2,bp1,bp2,inserts,max_ins):        
            span_bp1 = np.append(span_bp1,r1)
            span_bp2 = np.append(span_bp2,r2)
            rc['spanning'] = rc['spanning']+1
  
    if use_dir:
        rc['bp1_dir'] = bp1['dir']
        rc['bp2_dir'] = bp2['dir']
    else:
        # estimate SV directionality
        rc['bp1_dir'] = get_dir(split_bp1,span_bp1,loc1_reads,pos1,sc_len)
        rc['bp2_dir'] = get_dir(split_bp2,span_bp2,loc2_reads,pos2,sc_len)

    print('processed %d reads at loc1; %d reads at loc2' % (len(loc1_reads),len(loc2_reads)))
    return rc

def get_params(bam,mean_dp,max_cn,rlen,insert_mean,insert_std,out):
    inserts = [insert_mean,insert_std]
    if rlen<0:
        rlen = bamtools.estimateTagSize(bam)
    if inserts[0]<0 or inserts[1]<0:
        inserts = bamtools.estimateInsertSizeDistribution(bam)
    else:
        inserts[0] = inserts[0]+(rlen*2) #derive fragment size
    
    max_dp = ((mean_dp*(params.window*2))/rlen)*max_cn
    max_ins = 2*inserts[1]+inserts[0] #actually the *fragment* size

    with open('%s_params.txt'%out,'w') as outp:
        outp.write('read_len\tinsert_mean\tinsert_std\tinsert_cutoff\tmax_dep\n')
        outp.write('%d\t%f\t%f\t%f\t%d'%(rlen,inserts[0]-(rlen*2),inserts[1],max_ins-(rlen*2),max_dp))

    return rlen, inserts, max_dp, max_ins

def remove_duplicates(svs,use_dir):
    for idx,row in enumerate(svs):
        #reorder breakpoints based on position or chromosomes
        if use_dir:
            bp1_chr, bp1_pos, bp1_dir, bp2_chr, bp2_pos, bp2_dir = row
            if (bp1_chr!=bp2_chr and bp1_chr>bp2_chr) or (bp1_chr==bp2_chr and bp1_pos > bp2_pos):
                svs[idx] = (bp2_chr,bp2_pos,bp2_dir,bp1_chr,bp1_pos,bp1_dir)
        else:
            bp1_chr, bp1_pos, bp2_chr, bp2_pos = row
            if (bp1_chr!=bp2_chr and bp1_chr>bp2_chr) or (bp1_chr==bp2_chr and bp1_pos > bp2_pos):
                svs[idx] = (bp2_chr,bp2_pos,bp1_chr,bp1_pos)
    return np.unique(svs)

def load_input_vcf(svin):
    sv_dtype = [s for i,s in enumerate(params.sv_dtype) if i not in [2,5]]
    
    sv_vcf = vcf.Reader(filename=svin)
    sv_dict = OrderedDict()
    for sv in sv_vcf:
        
        if sv.FILTER is not None:
            if len(sv.FILTER)>0:
                continue
        
        sv_dict[sv.ID] = {'CHROM': sv.CHROM, 'POS': sv.POS, 'INFO': sv.INFO}

#    svs = OrderedDict()
#    sv_vcf = np.genfromtxt(svin,dtype=params.sv_vcf_dtype,delimiter='\t',comments="#")
#    keys = [key[0] for key in params.sv_vcf_dtype]
#    
#    for sv in sv_vcf:
#        sv_id = sv['ID']
#        svs[sv_id] = OrderedDict()
#        for key,sv_data in zip(keys,sv):
#            if key=='INFO' or key=='ID': continue
#            svs[sv_id][key] = sv_data
#         
#        info = map(methodcaller('split','='),sv['INFO'].split(';'))
#        svs[sv_id]['INFO'] = OrderedDict()
#        
#        for i in info:
#            if len(i)<2: continue
#            name = i[0]
#            data = i[1]
#            svs[sv_id]['INFO'][name] = data
    
    svs = np.empty(0,sv_dtype)
    procd = np.empty(0,dtype='S50')

    for sv_id in sv_dict:
        try:
            sv = sv_dict[sv_id]
            mate_id = sv['INFO']['MATEID']
            mate = sv_dict[mate_id]
            
            if (sv_id in procd) or (mate_id in procd): 
                continue
            
            bp1_chr = sv['CHROM']
            bp1_pos = sv['POS']
            bp2_chr = mate['CHROM']
            bp2_pos = mate['POS']
            
            procd = np.append(procd,[sv_id,mate_id])
            new_sv = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos)],dtype=sv_dtype)        
            svs = np.append(svs,new_sv)
        except KeyError:
            print("SV %s improperly paired or missing attributes"%sv_id)
            continue
    
    return svs

def load_input_socrates(svin,rlen,use_dir):
    sv_dtype =  [s for s in params.sv_dtype] if use_dir else [s for i,s in enumerate(params.sv_dtype) if i not in [2,5]]

    soc_in = np.genfromtxt(svin,delimiter='\t',names=True,dtype=None)
    svs = np.empty(0,dtype=sv_dtype)

    for row in soc_in:
        bp1 = row['C1_anchor'].split(':')
        bp2 = row['C1_realign'].split(':')
        bp1_chr, bp1_pos = bp1[0], int(bp1[1]) 
        bp2_chr, bp2_pos = bp2[0], int(bp2[1])
        #classification = row['classification']
        if not bp1_chr in params.valid_chroms or not bp2_chr in params.valid_chroms:
            continue
        if (bp1_chr==bp2_chr and abs(bp1_pos-bp2_pos)<(rlen*2)):
            continue
        if row['C1_avg_realign_mapq']<params.min_mapq or row['C2_avg_realign_mapq']<params.min_mapq:
            continue
        add_sv = np.empty(0)
        if use_dir:
            bp1_dir = row['C1_anchor_dir']
            bp2_dir = row['C1_realign_dir']
            add_sv = np.array([(bp1_chr,bp1_pos,bp1_dir,bp2_chr,bp2_pos,bp2_dir)],dtype=sv_dtype)
        else:
            add_sv = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos)],dtype=sv_dtype)
        svs = np.append(svs,add_sv)

    return remove_duplicates(svs,use_dir)

def load_input_simple(svin,use_dir):
    sv_dtype =  [s for s in params.sv_dtype] if use_dir else [s for i,s in enumerate(params.sv_dtype) if i not in [2,5]]

    sv_tmp = np.genfromtxt(svin,delimiter='\t',names=True,dtype=None)
    svs = np.empty(0,dtype=sv_dtype)
    for row in sv_tmp:
        bp1_chr = str(row['bp1_chr'])
        bp1_pos = int(row['bp1_pos'])
        bp2_chr = str(row['bp2_chr'])
        bp2_pos = int(row['bp2_pos'])
        add_sv = np.empty(0)
        if use_dir:
            bp1_dir = str(row['bp1_dir'])
            bp2_dir = str(row['bp2_dir'])
            add_sv = np.array([(bp1_chr,bp1_pos,bp1_dir,bp2_chr,bp2_pos,bp2_dir)],dtype=sv_dtype)
        else:
            add_sv = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos)],dtype=sv_dtype)
        svs = np.append(svs,add_sv)
    return remove_duplicates(svs,use_dir)

def proc_svs(args):

    svin         = args.svin
    bam          = args.bam
    out          = args.out
    mean_dp      = float(args.mean_dp)
    sc_len       = int(args.sc_len)
    max_cn       = int(args.max_cn)
    rlen         = int(args.rlen)
    insert_mean  = float(args.insert_mean)
    insert_std   = float(args.insert_std)
    simple       = args.simple_svs
    socrates     = args.socrates
    use_dir      = args.use_dir
   
    if not (simple or socrates): use_dir = False #vcfs don't have dirs

    db_out = '%s_svinfo.db'%out
    outf = '%s_svinfo.txt'%out
    
    dirname = os.path.dirname(out)
    if dirname!='' and not os.path.exists(dirname):
        os.makedirs(dirname)

    rlen, inserts, max_dp, max_ins = get_params(bam, mean_dp, max_cn, rlen, insert_mean, insert_std, out)

    # write header output
    header_out = ['ID'] + [h[0] for idx,h in enumerate(params.sv_dtype) if idx not in [2,5]] #don't include dir fields
    header_out.extend([h[0] for h in params.sv_out_dtype])
    
    with open('%s_svinfo.txt'%out,'w') as outf:        
        writer = csv.writer(outf,delimiter='\t',quoting=csv.QUOTE_NONE)
        writer.writerow(header_out)

    bp_dtype = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S1')] if use_dir \
                 else [('chrom','S20'),('start', int), ('end', int)]
    bp1_chr, bp1_pos, bp1_dir, bp2_chr, bp2_pos, bp2_dir = [h[0] for h in params.sv_dtype]
             
    svs = np.empty(0)
    if simple:
        svs = load_input_simple(svin,use_dir)
    elif socrates:
        svs = load_input_socrates(svin,rlen,use_dir)
    else:
        svs = load_input_vcf(svin)
    print("Extracting data from %d SVs"%len(svs))

    svd_prevResult, prevSV = None, None
    id = 0
    for row in svs:        
        sv_prop = row[bp1_chr],row[bp1_pos],row[bp2_chr],row[bp2_pos]
        sv_str = '%s:%d|%s:%d'%sv_prop

        print('processing %s'%sv_str)
        sv_rc = np.empty(0)
        if use_dir:
            bp1 = np.array((row[bp1_chr],row[bp1_pos]-params.window,row[bp1_pos]+params.window,row[bp1_dir]),dtype=bp_dtype)
            bp2 = np.array((row[bp2_chr],row[bp2_pos]-params.window,row[bp2_pos]+params.window,row[bp2_dir]),dtype=bp_dtype)
            sv_rc  = get_sv_read_counts(bp1,bp2,bam,inserts,max_dp,max_ins,sc_len,use_dir)
        else:
            bp1 = np.array((row[bp1_chr],row[bp1_pos]-params.window,row[bp1_pos]+params.window),dtype=bp_dtype)
            bp2 = np.array((row[bp2_chr],row[bp2_pos]-params.window,row[bp2_pos]+params.window),dtype=bp_dtype)
            sv_rc  = get_sv_read_counts(bp1,bp2,bam,inserts,max_dp,max_ins,sc_len,use_dir)
        
        if len(sv_rc) > 0:
            norm1 = int(sv_rc['bp1_split_norm']+sv_rc['bp1_span_norm'])
            norm2 = int(sv_rc['bp2_split_norm']+sv_rc['bp2_span_norm'])
            support = float(sv_rc['bp1_split'] + sv_rc['bp2_split'] + sv_rc['spanning'])
        
            sv_rc['norm1'] = norm1 
            sv_rc['norm2'] = norm2
            sv_rc['support'] = support
            sv_rc['vaf1'] = support / (support + norm1) if support!=0 else 0
            sv_rc['vaf2'] = support / (support + norm2) if support!=0 else 0
           
            #classify event
            sv = (row['bp1_chr'],row['bp1_pos'],sv_rc['bp1_dir'][0],row['bp2_chr'],row['bp2_pos'],sv_rc['bp2_dir'][0])
            svd_result = svd.detect(prevSV,svd_prevResult,sv)
            svd_prevResult,prevSV = svd_result,sv
            sv_rc['classification'] = svd.getResultType(svd_result)
            id = id if svd_result[0]==svd.SVtypes.interspersedInsertion else id+1
#            if row['classification']!=sv_rc['classification']:
#                with open('assign.txt','a') as outp:
#                    outp.write('%s:%d %s %s:%d %s'%sv + ' soc_class = %s; integrated = %s\n' % (row['classification'],sv_rc['classification']))            
            sv_out = [id] + [r for idx,r in enumerate(row) if idx not in [2,5]] if use_dir else \
                     [id] + [r for idx,r in enumerate(row)]
            sv_out.extend([rc for rc in sv_rc[0]])

            with open('%s_svinfo.txt'%out,'a') as outf:
                writer = csv.writer(outf,delimiter='\t',quoting=csv.QUOTE_NONE)
                writer.writerow(sv_out)
    #post-process: look for translocations
    sv_info = np.genfromtxt('%s_svinfo.txt'%out,delimiter='\t',names=True,dtype=None)
    trx_label    = svd.getResultType([svd.SVtypes.translocation])
    intins_label = svd.getResultType([svd.SVtypes.interspersedInsertion])
    rewrite=False
    for idx,sv in enumerate(sv_info):
        if sv['classification']==intins_label:
            rewrite=True
            sv_info[idx-1]['classification'] = intins_label
            translocs = svd.detectTransloc(idx,sv_info)
            if len(translocs)>0:
                for i in translocs: sv_info[i]['classification'] = trx_label
    if rewrite:
        with open('%s_svinfo.txt'%out,'w') as outf:
            writer = csv.writer(outf,delimiter='\t',quoting=csv.QUOTE_NONE)
            writer.writerow(header_out)
            for sv_out in sv_info:
                writer.writerow(sv_out)
