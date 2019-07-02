#!/usr/bin/env python
'''
Simple script to match output SVs across
two samples into one output file with
matching SVs on the same line and unmatched
SVs adjacent to blank records
'''

import argparse
import ipdb
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument('sv1')
parser.add_argument('sv2')
parser.add_argument('sample1')
parser.add_argument('sample2')
parser.add_argument('outfile')
args = parser.parse_args()

ff = 15 # fuzzy factor
sv1 = args.sv1
sv2 = args.sv2
outfile = args.outfile
sample1 = args.sample1
sample2 = args.sample2

def is_same_sv(sv1,sv2,threshold=20):
    sv1_chr1, sv1_bp1, sv1_dir1, sv1_chr2, sv1_bp2, sv1_dir2 = sv1
    sv2_chr1, sv2_bp1, sv2_dir1, sv2_chr2, sv2_bp2, sv2_dir2 = sv2

    if sv1_chr1==sv2_chr1 and sv1_chr2==sv2_chr2 and sv1_dir1 == sv2_dir1 and sv1_dir2 == sv2_dir2:
        if abs(sv1_bp1-sv2_bp1)<threshold and abs(sv1_bp2-sv2_bp2)<threshold:
            return True
    # same SVs with the order flipped
    if sv1_chr2==sv2_chr1 and sv1_chr1==sv2_chr2 and sv1_dir2 == sv2_dir1 and sv1_dir1 == sv2_dir2:
        if abs(sv1_bp2-sv2_bp1)<threshold and abs(sv1_bp1-sv2_bp2)<threshold:
            return True
    return False

def match_events(sv1,sv2):
    sv1_idxs = []
    sv2_idxs = []
    for idx1,s1 in sv1.iterrows():
        for idx2,s2 in sv2.iterrows():
            s1_tmp = [s1.bp1_chr,s1.bp1_pos,s1.bp1_dir,s1.bp2_chr,s1.bp2_pos,s1.bp2_dir]
            s2_tmp = [s2.bp1_chr,s2.bp1_pos,s2.bp1_dir,s2.bp2_chr,s2.bp2_pos,s2.bp2_dir]
            if is_same_sv(s1_tmp,s2_tmp,ff):
                sv1_idxs.append(idx1)
                sv2_idxs.append(idx2)
                #uset = np.append(uset,[[s1,s2]],axis=0)
                break

    return np.array(sv1_idxs),np.array(sv2_idxs)

if __name__ == "__main__":
    sv1_info = pd.DataFrame(pd.read_csv(sv1,delimiter='\t',dtype=None, low_memory=False))
    sv2_info = pd.DataFrame(pd.read_csv(sv2,delimiter='\t',dtype=None, low_memory=False))
    sv1_info.bp1_chr = sv1_info.bp1_chr.map(str)
    sv1_info.bp2_chr = sv1_info.bp2_chr.map(str)
    sv2_info.bp1_chr = sv2_info.bp1_chr.map(str)
    sv2_info.bp2_chr = sv2_info.bp2_chr.map(str)

    sv1_idxs,sv2_idxs = match_events(sv1_info,sv2_info)
    sv1_info.columns = ['%s_%s' % (sample1,col) for col in sv1_info.columns.values]
    sv2_info.columns = ['%s_%s' % (sample2,col) for col in sv2_info.columns.values]
    sv1_match = sv1_info.loc[sv1_idxs]
    sv2_match = sv2_info.loc[sv2_idxs]

    #fields = ['ID','bp1_chr','bp1_pos','bp1_dir','bp2_chr','bp2_pos','bp2_dir','adjusted_vaf','classification']
    #sv_match = sv1_match[fields]

    sv_match = pd.concat([sv1_match.reset_index(), sv2_match.reset_index()], axis=1)
    del sv_match['index']

    sv1_unmatch_idx = [idx for idx in sv1_info.index.values if idx not in sv1_idxs]
    sv2_unmatch_idx = [idx for idx in sv2_info.index.values if idx not in sv2_idxs]
    sv1_unmatch = sv1_info.loc[sv1_unmatch_idx]
    sv2_unmatch = sv2_info.loc[sv2_unmatch_idx]

    sv1_unmatch = sv1_unmatch.join(pd.DataFrame(dtype=sv2_match.dtypes, columns=sv2_match.columns))
    sv2_unmatch = sv2_unmatch.join(pd.DataFrame(dtype=sv1_match.dtypes, columns=sv1_match.columns))
    sv2_unmatch = sv2_unmatch[sv1_unmatch.columns.values]
    sv_match = sv_match.append(sv1_unmatch)
    sv_match = sv_match.append(sv2_unmatch)

    sv_match.to_csv(outfile,sep='\t',index=False,na_rep='')



