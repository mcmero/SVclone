'''
Takes an SV matrix and attempts to build a phylogeny
'''
import subprocess
import ipdb
import sys
import os
import shutil
import pandas as pd

from . import cluster
from . import parameters as param
# from ete2 import Tree

import numpy as np
from numpy import loadtxt
from numpy import zeros
from scipy import stats

def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def merge_clusters(df,clus_info,clus_merged,clus_members,cparams):
    if len(clus_info)==1:
        return (clus_info,clus_members)

    to_del = []
    for idx,ci in clus_info.iterrows():
        if idx+1 >= len(clus_info):
            clus_merged.loc[idx] = clus_info.loc[idx]
            break
        cn = clus_info.loc[idx+1]
        if abs(ci.phi - float(cn.phi)) < param.subclone_diff:
            print("\nReclustering similar clusters...")
            new_size = ci['size'] + cn['size']
            new_members = np.concatenate([clus_members[idx],clus_members[idx+1]])
            new_phi = cluster.recluster(df.loc[new_members], cparams[0], cparams[1], cparams[2])#, 7000)
            clus_merged.loc[idx] = np.array([new_size,ci.clus_id,new_phi])
            clus_members[idx] = new_members
            to_del.append(idx+1)
            print('\n')
            if idx+2 < len(clus_info):
                clus_merged.loc[idx+2:] = clus_info.loc[idx+2:]
            break
        else:
            clus_merged.loc[idx] = clus_info.loc[idx]
    
    if len(to_del)>0:
        clus_members = np.delete(clus_members,to_del)
        # clean and reorder merged dataframe
        clus_merged = clus_merged[pd.notnull(clus_merged['clus_id'])]
        clus_merged = clus_merged.sort('phi',ascending=False)
        clus_merged.index = range(len(clus_merged))
        print("\nMerged clusters")
        print('Input')
        print(clus_info)
        print('Merged')
        print(clus_merged)

    if len(clus_info) == len(clus_merged):
        return (clus_merged,clus_members)
    else: 
        new_df = pd.DataFrame(columns=clus_merged.columns,index=clus_merged.index)
        return merge_clusters(df,clus_merged,new_df,clus_members,cparams)


def infer_tree(df,samples,pi,rlen,insert,out):
    mcmc = cluster.cluster(df,pi,rlen,insert)
    npoints = len(df.spanning.values)

    # assign points to highest probabiliy cluster
    z_trace = mcmc.trace('z')[:]
    clus_counts = [np.bincount(z_trace[:,i]) for i in range(npoints)]
    clus_max_prob = [index_max(c) for c in clus_counts]
    clus_mp_counts = np.bincount(clus_max_prob)
    clus_idx = np.nonzero(clus_mp_counts)[0]
    clus_mp_counts = clus_mp_counts[clus_idx]
    
    # cluster distribution
    clus_info = pd.DataFrame(clus_idx,columns=['clus_id'])
    clus_info['size'] = clus_mp_counts

    # get cluster means
    center_trace = mcmc.trace("phi_k")[:]
    # param.plot_clusters(center_trace,npoints,clus_idx)
    
    phis = np.array([np.mean(center_trace[:,cid]) for cid in clus_info.clus_id.values])
    clus_info['phi'] = phis
    clus_info = clus_info.sort('phi',ascending=False)

    # filter clones below thresholds
    clus_info = clus_info[clus_info['phi'].values>param.subclone_threshold]
    clus_info = clus_info[sum(clus_info['size'].values)*param.subclone_sv_prop<clus_info['size'].values]
    clus_info.index = range(len(clus_info))

    if len(clus_info) < 1:
        print("Could not converge on any major clusters! Exiting.")
        return None
   
    cparams = [pi,rlen,insert] #params required for clustering
    clus_members = np.array([np.where(np.array(clus_max_prob)==i)[0] for i in clus_info.clus_id.values])    
    clus_init = pd.DataFrame(columns=clus_info.columns,index=clus_info.index)
    clus_new,clus_members = merge_clusters(df,clus_info,clus_init,clus_members,cparams)

    if (clus_info.phi.values[0]*2.)>(1+param.subclone_diff):
        print("\nWarning! VAF cluster exceeds 0.5. Germline breaks may not have been adequately filtered " + \
              "or provided copy-numbers are inaccurate. Tumour purity may also be over-estimated. " +
              "Please double-check and rerun.")
    
    print("\nFiltered & merged clusters")
    clus_new.clus_id = map(int,clus_new.clus_id.values)
    clus_new['size'] = map(int,clus_new['size'].values)
    print(clus_new)

    clus_new.to_csv('%s.txt'%out,sep='\t',index=False)
    for idx,clus in clus_new.iterrows():
        cm = clus_members[idx]
        clus_df = df.loc[cm]
        clus_df.to_csv('%s_clus%d.txt'%(out,clus.clus_id),sep='\t',index=False)

