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
sc_tr = 25 #soft-clip threshold by which we call split reads
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
       
def get_normal_read_counts(chrom, start, end, bdir, bamf):
    outfile = pysam.AlignmentFile('sclipped.bam', "wh", header=bamf.header)
    pos = (start + end) / 2
    rc = pd.Series(np.zeros(4),index=['split_norm','span_norm','sc_reads','anomalous'])
    reg = ':'.join([str(chrom), str(start), str(end)])
    iter = bamf.fetch(region=reg)
    for x in iter:
        m = bamf.mate(x)
        #if is_soft_clipped(x) or is_soft_clipped(bamf.mate(x)):
        #    print (x.query_name)
        #if x.query_name=='SimSeq_30cc4c':
        #    ipdb.set_trace()
        if is_aligned_across_break(x,pos):
            rc.split_norm = rc.split_norm+1
            outfile.write(x)
            #print ('%s aligns across break' % x.reference_start)
        elif spans_across_break(x,bamf.mate(x),pos):
            rc.span_norm = rc.span_norm+1
            #print ('%s spans across break' % x.reference_start)
        elif is_soft_clipped(x):
            rc.sc_reads = rc.sc_reads+1
    #rc.span_norm = rc.span_norm/2
    print rc
    outfile.close()

def run(svin,bam):    
    bamf = pysam.AlignmentFile(bam, "rb")
    svs = pd.read_csv(svin,delimiter='\t')
    for idx,row in svs.iterrows():
        bp1_start = row.bp1_pos-window
        bp1_end = row.bp1_pos+window
        print row.bp1_pos
        get_normal_read_counts(row.bp1_chr,bp1_start,bp1_end,row.bp1_dir,bamf)
        break
    bamf.close()

