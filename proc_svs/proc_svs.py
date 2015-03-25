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
    if r.query_alignment_start == 0:
        return r.reference_end > (pos - tr) and r.reference_end < (pos + tr) and \
            (r.query_length - r.query_alignment_end)
    else:
        return r.reference_start > (pos - tr) and r.reference_start < (pos + tr) and \
            (r.query_alignment_start >= sc_len)

#def is_supporting_spanning_pair(r,m,pos):
     

def get_read_counts(chrom, start, end, bamf, outfile):
    pos = (start + end) / 2
    rc = pd.Series(np.zeros(5),index=['split_norm','span_norm','split','spanning','anomalous'])
    reg = ':'.join([str(chrom), str(start), str(end)])
    iter = bamf.fetch(region=reg)
    for x in iter:
        #TODO: check whether read has a congruent insert size
        m = bamf.mate(x)
        if is_aligned_across_break(x,pos):
            rc.split_norm = rc.split_norm+1
            #print ('%s aligns across break' % x.reference_start)
        elif spans_across_break(x,m,pos):
            rc.span_norm = rc.span_norm+1
            #print ('%s spans across break' % x.reference_start)
        elif is_supporting_split_read(x,pos):
            rc.split = rc.split+1
        else:
            outfile.write(x)
        #if x.query_name=='SimSeq_ee41f':
        #    is_supporting_split_read(x,pos)
    return rc

def run(svin,bam):    
    bamf = pysam.AlignmentFile(bam, "rb")
    svs = pd.read_csv(svin,delimiter='\t')
    columns = ['bp1_split_norm','bp1_span_norm','bp1_split','bp1_spanning','bp1_anomalous', \
               'bp2_split_norm','bp2_span_norm','bp2_split','bp2_spanning','bp2_anomalous']
    rinfo = pd.DataFrame(columns=columns,index=range(len(svs.index)))
    for idx,row in svs.iterrows():
        bp1_start = row.bp1_pos-window
        bp1_end = row.bp1_pos+window
        bp2_start = row.bp2_pos-window
        bp2_end = row.bp2_pos+window
        
        #print row.bp1_pos
        outfile = pysam.AlignmentFile('anom.bam', "wh", header=bamf.header)
        bp1 = get_read_counts(row.bp1_chr,bp1_start,bp1_end,bamf,outfile)
        bp2 = get_read_counts(row.bp2_chr,bp2_start,bp2_end,bamf,outfile)
        
        for i in range(len(bp1)):
            bp1.index.values[i] = 'bp1_%s' % bp1.index.values[i]
            bp2.index.values[i] = 'bp2_%s' % bp2.index.values[i]
        #ipdb.set_trace()
        rinfo.loc[idx] = pd.concat([bp1,bp2])
        print rinfo
    print rinfo
    bamf.close()

