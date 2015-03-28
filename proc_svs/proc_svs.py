'''
Get SV output, run post processing
'''

import os
import csv
import string
import numpy as np
import subprocess
import argparse
import itertools
import ipdb
import pysam
import pandas as pd
import bamtools
from subprocess import call
from itertools import permutations
from numpy import loadtxt
from scipy import stats

tr = 5 #threshold by how much read has to overlap breakpoint
sc_len = 25 #soft-clip threshold by which we call split reads
window = 300

read_dtype = [('query_name', 'S25'), ('chrom', str), ('ref_start', int), ('ref_end', int), \
              ('align_start', int), ('align_end', int), ('len', int), ('is_reverse', np.bool)]

def read_to_array(x,bamf):
    chrom = bamf.getrname(x.reference_id)
    read = np.array((x.query_name,chrom,x.reference_start,x.reference_end,x.query_alignment_start,
                     x.query_alignment_end,x.query_length,np.bool(x.is_reverse)),dtype=read_dtype)
    return read

#def is_aligned_across_break(r,pos):
#   return r.reference_start < (pos - tr) and r.reference_end > (pos + tr) and \
#       r.query_alignment_start < (tr*2) and (r.query_alignment_end + r.reference_start - r.reference_end < (tr*2))


#def is_soft_clipped(r):
#    return r.query_alignment_start != 0 or (r.query_length + r.reference_start != r.reference_end)

def is_soft_clipped(r):
    return r['align_start'] != 0 or (r['len'] + r['ref_start'] != r['ref_end'])

def is_normal_across_break(r,pos):
   return r['ref_start'] < (pos - tr) and r['ref_end'] > (pos + tr) and \
       r['align_start'] < (tr*2) and (r['align_end'] + r['ref_start'] - r['ref_end'] < (tr*2))

def is_normal_spanning_break(r,m,pos):
    if not (is_soft_clipped(r) or is_soft_clipped(m)):
        if (not r['is_reverse'] and m['is_reverse']) or (r['is_reverse'] and not m['is_reverse']):
            return r['ref_end'] < (pos + tr) and m['ref_start'] > (pos - tr)
    return False

#def is_supporting_split_read(r,pos):
#    """
#    return whether read is a supporting split read
#    doesn't yet check whether the soft-clip aligns
#    to the other side - should it, or is that implicit?
#    """
#    if r.query_alignment_start == 0:
#        return r.reference_end > (pos - tr) and r.reference_end < (pos + tr) and \
#            (r.query_length - r.query_alignment_end)
#    else:
#        return r.reference_start > (pos - tr) and r.reference_start < (pos + tr) and \
#            (r.query_alignment_start >= sc_len)

def is_supporting_split_read(r,pos):
    """
    return whether read is a supporting split read
    doesn't yet check whether the soft-clip aligns
    to the other side - should it, or is that implicit?
    """
    if r['align_start'] == 0:
        return r['ref_end'] > (pos - tr) and r['ref_end'] < (pos + tr) and \
            (r['len'] - r['align_end'])
    else:
        return r['ref_start'] > (pos - tr) and r['ref_start'] < (pos + tr) and \
            (r['align_start'] >= sc_len)

def is_supporting_spanning_pair(r,m,bp1,bp2,sv_type,inserts):
    pos1 = (bp1['start'] + bp1['end']) / 2
    pos2 = (bp2['start'] + bp2['end']) / 2
    dir1 = bp1['dir']
    dir2 = bp2['dir']

    #check if pair is facing the right way (only if it is a duplication)
    if sv_type == 'DUP':
        #if dup is the "correct way around" i.e. pos1 < pos2
        if dir1=='-' and dir2=='+':
            if not (not r['is_reverse'] and m['is_reverse']):
                return False
        #dup is wrong way round
        elif dir1=='+' and dir2=='-': 
            if not (r['is_reverse'] and not m['is_reverse']):
                return False
        elif dir1 == dir2:
            return False #not a duplicaton?
        
    ins_dist = -1    
    if dir1=='-':
        ins_dist = pos1 - r['ref_start']
    else:
        ins_dist = r['ref_end'] - pos1

    if dir2=='-':
        ins_dist = ins_dist + (pos1 - r['ref_start'])
    else:
        ins_dist = ins_dist + (r['ref_end'] - pos1)

    if ins_dist>0:
        #print m
        #print r
        return ins_dist < (2*inserts[1]+inserts[0])
    else:
        return False

def get_loc_reads(bp,bamf):
    loc = '%s:%d:%d' % (bp['chrom'], bp['start'], bp['end'])
    
    loc_reads = np.empty([0,len(read_dtype)],dtype=read_dtype)    
    iter_loc = bamf.fetch(region=loc)
    for x in iter_loc:
        read = read_to_array(x,bamf) 
        loc_reads = np.append(loc_reads,read)
    
    loc_reads = np.sort(loc_reads,axis=0,order=['query_name','ref_start'])
    return loc_reads


def get_sv_read_counts(bp1,bp2,bamf,columns,sv_type,inserts):
    #ipdb.set_trace()
    pos1 = (bp1['start'] + bp1['end']) / 2
    pos2 = (bp2['start'] + bp2['end']) / 2
    loc1_reads = get_loc_reads(bp1,bamf)
    loc2_reads = get_loc_reads(bp2,bamf)
    
    rc = pd.Series(np.zeros(len(columns)),index=columns)
    reproc = np.empty([0,len(read_dtype)],dtype=read_dtype)
    
    for idx,x in enumerate(loc1_reads):
        r1 = loc1_reads[idx]
        r2 = loc1_reads[idx+1] if (idx+2)<=len(loc1_reads) else None
        if is_normal_across_break(x,pos1):
            rc['bp1_split_norm'] = rc['bp1_split_norm']+1 
        elif is_supporting_split_read(x,pos1):
            rc['bp1_split'] = rc['bp1_split']+1 
        elif r2!=None and r1['query_name']==r2['query_name'] and is_normal_spanning_break(r1,r2,pos1):
                rc['bp1_span_norm'] = rc['bp1_span_norm']+1 
        else:
            reproc = np.append(reproc,x) #may be spanning support or anomalous
    
    for idx,x in enumerate(loc2_reads):
        r1 = loc2_reads[idx]
        r2 = loc2_reads[idx+1] if (idx+2)<=len(loc2_reads) else None
        if is_normal_across_break(x,pos2):
            rc['bp2_split_norm'] = rc['bp2_split_norm']+1 
        elif is_supporting_split_read(x,pos2):
            rc['bp2_split'] = rc['bp2_split']+1 
        elif r2!=None and r1['query_name']==r2['query_name'] and is_normal_spanning_break(r1,r2,pos2):
                rc['bp2_span_norm'] = rc['bp2_span_norm']+1 
        else:
            reproc = np.append(reproc,x) #may be spanning support or anomalous
     
    reproc = np.sort(reproc,axis=0,order=['query_name','ref_start'])
    for idx,x in enumerate(reproc):
        if idx+1 >= len(reproc):
            break        
        if reproc[idx+1]['query_name'] != reproc[idx]['query_name']:
            #not paired
            rc['bp1_anomalous'] = rc['bp1_anomalous']+1 
            continue
        mate = reproc[idx+1]
        r1 = x
        r2 = mate
        if bp1['chrom']!=bp2['chrom'] and x['chrom']==bp2['chrom']:
            r1 = mate
            r2 = x
        elif pos1 > pos2 and bp1['chrom']==bp2['chrom']:
            r1 = mate
            r2 = x
        if is_supporting_spanning_pair(r1,r2,bp1,bp2,sv_type,inserts):
            rc['spanning'] = rc['spanning']+1 
        else:
            rc['bp1_anomalous'] = rc['bp1_anomalous']+1
    
    return rc

#def get_read_counts(chrom, start, end, bamf):
#    pos = (start + end) / 2
#    anom = np.empty([0,len(read_dtype)],dtype=read_dtype)    
#    rc = pd.Series(np.zeros(4),index=['split_norm','span_norm','split','anomalous'])
#    reg = ':'.join([str(chrom), str(start), str(end)])
#    iter = bamf.fetch(region=reg)
#    for x in iter:
#        m = bamf.mate(x)
#        #ipdb.set_trace()
#        if is_aligned_across_break(x,pos):
#            rc.spl'it_norm = rc.split_norm+1
#            #print ('%s aligns across break' % x.reference_start)
#        elif spans_across_break(x,m,pos):
#            rc.span_norm = rc.span_norm+1
#            #print ('%s spans across break' % x.reference_start)
#        elif is_supporting_split_read(x,pos):
#            rc.split = rc.split+1
#        else:
#            rc.anomalous = rc.anomalous+1
#            read = read_to_array(x,bamf) 
#            anom = np.append(anom,read) 
#            #tmp_bam.write(x)
#
#        #if x.query_name=='SimSeq_ee41f':
#        #    is_supporting_split_read(x,pos)
#    return rc, anom

#def get_spanning_support(bp1_chr, bp1_start, bp1_end, bp1_dir, bp2_chr, bp2_start, bp2_end, bp2_dir, sv_type, anom, inserts):
#    #tmp = pysam.AlignmentFile(tmp_bam, "r")
#    pos1 = (bp1_start + bp1_end) / 2
#    pos2 = (bp2_start + bp2_end) / 2
#    span = 0
#    
#    #dtype = [('query_name', 'S25'), ('ref_start', int), ('ref_end', int), \
#    #         ('align_start', int), ('align_end', int), ('is_reverse', np.bool)]
#    #span_arr = np.empty([0,6],dtype=dtype)    
#    #iter = tmp.fetch()
#    #for x in iter:
#    #    if is_soft_clipped(x):
#    #        continue
#    #    read = np.array((x.query_name,x.reference_start,x.reference_end, 
#    #                    x.query_alignment_start,x.query_alignment_end,
#    #                   np.bool(x.is_reverse)),dtype=dtype)
#    #    span_arr = np.append(span_arr,read) 
#    #span_arr = np.sort(span_arr,axis=0,order='query_name')
#    
#    for idx,r in enumerate(anom):
#        if idx+1 >= len(anom):
#            break
#        if anom[idx+1]['query_name'] != anom[idx]['query_name']:
#            #could not find mate
#            continue
#
#        #code here that determines which read is r1 (paired to bp1) and which to r2        
#        #TODO - must be a better way to do this
#        if bp1_chr!=bp2_chr and anom[idx]['chrom']==bp1_chr:
#            if is_supporting_spanning_pair(anom[idx],anom[idx+1],bp1_chr,bp2_chr,pos1,pos2,bp1_dir,bp2_dir,sv_type,inserts):
#                span = span + 1
#        elif bp1_chr!=bp2_chr and anom[idx]['chrom']==bp2_chr:
#            if is_supporting_spanning_pair(anom[idx+1],anom[idx],bp1_chr,bp2_chr,pos1,pos2,bp1_dir,bp2_dir,sv_type,inserts):
#                span = span + 1
#        elif bp1_chr==bp2_chr:            
#            if pos1 > pos2:
#                if is_supporting_spanning_pair(anom[idx+1],anom[idx],bp1_chr,bp2_chr,pos1,pos2,bp1_dir,bp2_dir,sv_type,inserts):
#                    span = span + 1
#            elif is_supporting_spanning_pair(anom[idx],anom[idx+1],bp1_chr,bp2_chr,pos1,pos2,bp1_dir,bp2_dir,sv_type,inserts):
#                span = span + 1
#    return span
#

def run(svin,bam,out):    
    inserts = bamtools.estimateInsertSizeDistribution(bam)
    
    bp_dtype = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S2')]

    bamf = pysam.AlignmentFile(bam, "rb")
    svs = pd.read_csv(svin,delimiter='\t')

    columns = ['bp1_split_norm','bp1_span_norm','bp1_split','bp1_anomalous', \
               'bp2_split_norm','bp2_span_norm','bp2_split','bp2_anomalous','spanning']
    rinfo = pd.DataFrame(columns=columns,index=range(len(svs.index)))
    
    for idx,row in svs.iterrows():        
        bp1 = np.array((row.bp1_chr,row.bp1_pos-window,row.bp1_pos+window,row.bp1_dir),dtype=bp_dtype)
        bp2 = np.array((row.bp2_chr,row.bp2_pos-window,row.bp2_pos+window,row.bp2_dir),dtype=bp_dtype)                     
        
        #print row.bp1_pos
        #tmp_bam = pysam.AlignmentFile('tmp.bam', "wh", header=bamf.header)
        sv_rc  = get_sv_read_counts(bp1,bp2,bamf,columns,row.classification,inserts)
        #bp2  = get_sv_read_counts(row.bp2_chr,bp2_start,bp2_end,bamf)
        #tmp_bam.close()

        #anom = np.append(anom1,anom2)
        #anom = np.sort(anom,axis=0,order=['query_name','ref_start'])

        #for i in range(len(bp1)):
        #    bp1.index.values[i] = 'bp1_%s' % bp1.index.values[i]
        #    bp2.index.values[i] = 'bp2_%s' % bp2.index.values[i]
        #ipdb.set_trace(i)
        #span = get_spanning_support(row.bp1_chr,bp1_start,bp1_end,bp1_dir,
                                    #row.bp2_chr,bp2_start,bp2_end,bp2_dir,
                                    #row.classification,anom,inserts)

        #newr = pd.concat([bp1,bp2])
        #newr['spanning'] = span
        rinfo.loc[idx] = sv_rc
        #print rinfo
        #break
    bamf.close()    
    rinfo = svs.join(rinfo)
    print rinfo
    rinfo.to_csv(out,sep="\t")    

