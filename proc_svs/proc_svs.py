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

def is_aligned_across_break(r,pos):
   return r.reference_start < (pos - tr) and r.reference_end > (pos + tr) and \
       r.query_alignment_start < (tr*2) and (r.query_alignment_end + r.reference_start - r.reference_end < (tr*2))

def is_soft_clipped(r):
    return r.query_alignment_start != 0 or (r.query_length + r.reference_start != r.reference_end)

def spans_across_break(r,m,pos):
    if not (is_soft_clipped(r) or is_soft_clipped(m)):
        if (not r.is_reverse and r.mate_is_reverse) or (r.is_reverse and not r.mate_is_reverse):
            return r.reference_end < (pos + tr) and m.reference_start > (pos - tr)
    return False

def is_supporting_split_read(r,pos):
    """
    return whether read is a supporting split read
    doesn't yet check whether the soft-clip aligns
    to the other side - should it, or is that implicit?
    """
    if r.query_alignment_start == 0:
        return r.reference_end > (pos - tr) and r.reference_end < (pos + tr) and \
            (r.query_length - r.query_alignment_end)
    else:
        return r.reference_start > (pos - tr) and r.reference_start < (pos + tr) and \
            (r.query_alignment_start >= sc_len)

def is_supporting_spanning_pair(r,m,chr1,chr2,pos1,pos2,dir1,dir2,sv_type,inserts):
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
        
    #TODO: quick dirty hack, fix later, properly
    #ipdb.set_trace() 

    dist1 = min([abs(pos1-r['ref_start']),abs(pos1-m['ref_start']),abs(pos1-r['ref_end']),abs(pos1-m['ref_end'])])
    dist2 = min([abs(pos2-r['ref_start']),abs(pos2-m['ref_start']),abs(pos2-r['ref_end']),abs(pos2-m['ref_end'])])
    if (dist1 + dist2) < (2*inserts[1]+inserts[0]):
        print m
        print r
        return True
    else:
        return False

def get_read_counts(chrom, start, end, bamf):
    pos = (start + end) / 2
    dtype = [('query_name', 'S25'), ('chrom', str), ('ref_start', int), ('ref_end', int), \
             ('align_start', int), ('align_end', int), ('is_reverse', np.bool)]
    anom = np.empty([0,7],dtype=dtype)    
    rc = pd.Series(np.zeros(4),index=['split_norm','span_norm','split','anomalous'])
    reg = ':'.join([str(chrom), str(start), str(end)])
    iter = bamf.fetch(region=reg)
    for x in iter:
        m = bamf.mate(x)
        #ipdb.set_trace()
        if is_aligned_across_break(x,pos):
            rc.split_norm = rc.split_norm+1
            #print ('%s aligns across break' % x.reference_start)
        elif spans_across_break(x,m,pos):
            rc.span_norm = rc.span_norm+1
            #print ('%s spans across break' % x.reference_start)
        elif is_supporting_split_read(x,pos):
            rc.split = rc.split+1
        else:
            rc.anomalous = rc.anomalous+1
            chrom = bamf.getrname(x.reference_id)
            read = np.array((x.query_name,chrom,x.reference_start,
                             x.reference_end,x.query_alignment_start,
                             x.query_alignment_end,np.bool(x.is_reverse)),dtype=dtype)
            anom = np.append(anom,read) 
            #tmp_bam.write(x)

        #if x.query_name=='SimSeq_ee41f':
        #    is_supporting_split_read(x,pos)
    return rc, anom

def get_spanning_support(bp1_chr, bp1_start, bp1_end, bp1_dir, bp2_chr, bp2_start, bp2_end, bp2_dir, sv_type, anom, inserts):
    #tmp = pysam.AlignmentFile(tmp_bam, "r")
    pos1 = (bp1_start + bp1_end) / 2
    pos2 = (bp2_start + bp2_end) / 2
    span = 0
    
    #dtype = [('query_name', 'S25'), ('ref_start', int), ('ref_end', int), \
    #         ('align_start', int), ('align_end', int), ('is_reverse', np.bool)]
    #span_arr = np.empty([0,6],dtype=dtype)    
    #iter = tmp.fetch()
    #for x in iter:
    #    if is_soft_clipped(x):
    #        continue
    #    read = np.array((x.query_name,x.reference_start,x.reference_end, 
    #                    x.query_alignment_start,x.query_alignment_end,
    #                   np.bool(x.is_reverse)),dtype=dtype)
    #    span_arr = np.append(span_arr,read) 
    #span_arr = np.sort(span_arr,axis=0,order='query_name')
    
    for idx,r in enumerate(anom):
        if idx+1 >= len(anom):
            break
        if anom[idx+1]['query_name'] != anom[idx]['query_name']:
            #could not find mate
            continue

        #code here that determines which read is r1 (paired to bp1) and which to r2        
        #TODO - must be a better way to do this
        if bp1_chr!=bp2_chr and anom[idx]['chrom']==bp1_chr:
            if is_supporting_spanning_pair(anom[idx],anom[idx+1],bp1_chr,bp2_chr,pos1,pos2,bp1_dir,bp2_dir,sv_type,inserts):
                span = span + 1
        elif bp1_chr!=bp2_chr and anom[idx]['chrom']==bp2_chr:
            if is_supporting_spanning_pair(anom[idx+1],anom[idx],bp1_chr,bp2_chr,pos1,pos2,bp1_dir,bp2_dir,sv_type,inserts):
                span = span + 1
            
        if is_supporting_spanning_pair(anom[idx],anom[idx+1],bp1_chr,bp2_chr,pos1,pos2,bp1_dir,bp2_dir,sv_type,inserts):
            span = span + 1
    return span


def run(svin,bam,out):    
    inserts = bamtools.estimateInsertSizeDistribution(bam)
    
    bamf = pysam.AlignmentFile(bam, "rb")
    svs = pd.read_csv(svin,delimiter='\t')
    columns = ['bp1_split_norm','bp1_span_norm','bp1_split','bp1_anomalous', \
               'bp2_split_norm','bp2_span_norm','bp2_split','bp2_anomalous','spanning']
    rinfo = pd.DataFrame(columns=columns,index=range(len(svs.index)))
    
    for idx,row in svs.iterrows():        
        bp1_start = row.bp1_pos-window
        bp1_end = row.bp1_pos+window
        bp1_dir = row.bp1_dir
        bp2_start = row.bp2_pos-window
        bp2_end = row.bp2_pos+window
        bp2_dir = row.bp2_dir
        
        #print row.bp1_pos
        #tmp_bam = pysam.AlignmentFile('tmp.bam', "wh", header=bamf.header)
        bp1,anom1  = get_read_counts(row.bp1_chr,bp1_start,bp1_end,bamf)
        bp2,anom2  = get_read_counts(row.bp2_chr,bp2_start,bp2_end,bamf)
        #tmp_bam.close()

        anom = np.append(anom1,anom2)
        anom = np.sort(anom,axis=0,order=['query_name','ref_start'])

        for i in range(len(bp1)):
            bp1.index.values[i] = 'bp1_%s' % bp1.index.values[i]
            bp2.index.values[i] = 'bp2_%s' % bp2.index.values[i]
        #ipdb.set_trace(i)
        span = get_spanning_support(row.bp1_chr,bp1_start,bp1_end,bp1_dir,
                                    row.bp2_chr,bp2_start,bp2_end,bp2_dir,
                                    row.classification,anom,inserts)

        newr = pd.concat([bp1,bp2])
        newr['spanning'] = span
        rinfo.loc[idx] = newr
        print rinfo
        break
    bamf.close()
    rinfo = svs.join(rinfo)
    rinfo.to_csv(out,sep="\t")    

