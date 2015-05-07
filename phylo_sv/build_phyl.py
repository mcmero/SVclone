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

def merge_clusters(df,clus_info,clus_merged,clus_members,cparams,out):
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
            new_phi,trace = cluster.recluster(df.loc[new_members], cparams[0], cparams[1], cparams[2])
            clus_members[idx] = new_members
            to_del.append(idx+1)
            
            df_trace = pd.DataFrame(np.transpose(trace[:]))
            df_trace.to_csv('%s/reclus%s_phi_trace.txt'%(out,ci.clus_id),sep='\t',index=False)
            clus_merged.loc[idx] = np.array([ci.clus_id,new_size,new_phi])
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
        return merge_clusters(df,clus_merged,new_df,clus_members,cparams,out)

def run_clust(df,pi,rlen,insert):
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
    center_trace = center_trace[len(center_trace)/4:]
    # param.plot_clusters(center_trace,npoints,clus_idx)
    
    phis = np.array([np.mean(center_trace[:,cid]) for cid in clus_info.clus_id.values])
    clus_info['phi'] = phis
    clus_info = clus_info.sort('phi',ascending=False)

    # filter clones below thresholds
    clus_info = clus_info[clus_info['phi'].values>param.subclone_threshold]
    clus_info = clus_info[sum(clus_info['size'].values)*param.subclone_sv_prop<clus_info['size'].values]
    clus_info.index = range(len(clus_info))
    
    clus_members = np.array([np.where(np.array(clus_max_prob)==i)[0] for i in clus_info.clus_id.values])    
    return clus_info,center_trace,z_trace,clus_members

def dump_trace(clus_info,center_trace,outf):
    traces = np.array([center_trace[:,cid] for cid in clus_info.clus_id.values])
    df_traces = pd.DataFrame(np.transpose(traces),columns=clus_info.clus_id)
    df_traces.to_csv(outf,sep='\t',index=False)

def infer_subclones(df,pi,rlen,insert,out,maxiters):

    clus_info,center_trace,ztrace,clus_members = None,None,None,None

    for i in range(maxiters):
        print("Cluster run: %d"%i)
        clus_info,center_trace,z_trace,clus_members = run_clust(df,pi,rlen,insert)
        sv_loss = 1-(sum(clus_info['size'])/float(len(df)))

        if len(clus_info) < 1:
            print("\nWarning! Could not converge on any major clusters.")
            continue
        elif sv_loss > param.tolerate_svloss:
            print("\nWarning! Lost %f of SVs." % sv_loss)
        elif (clus_info.phi.values[0]*2.)>(1+param.subclone_diff):
            print("\nWarning! Largest VAF cluster exceeds 0.5.")
   
        clus_out_dir = '%s/run%d'%(out,i)
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)
        
        dump_trace(clus_info,center_trace,'%s/phi_trace.txt'%clus_out_dir)
        dump_trace(clus_info,z_trace,'%s/z_trace.txt'%clus_out_dir)
        clus_info.to_csv('%s/clusters_premerged.txt'%(clus_out_dir),sep='\t',index=False)

        clus_init = pd.DataFrame(columns=clus_info.columns,index=clus_info.index)
        clus_info,clus_new_members = \
            merge_clusters(df,clus_info,clus_init,clus_members.copy(),[pi,rlen,insert],clus_out_dir)
        
        print("\nFiltered & merged clusters")
        clus_info.clus_id = map(int,clus_info.clus_id.values)
        clus_info['size'] = map(int,clus_info['size'].values)
        print(clus_info)

        clus_info.to_csv('%s/clusters.txt'%(clus_out_dir),sep='\t',index=False)
        for idx,clus in clus_info.iterrows():
            cm = clus_new_members[idx]
            clus_df = df.loc[cm]
            clus_df.to_csv('%s/clus%d_svs.txt'%(clus_out_dir,clus.clus_id),sep='\t',index=False)

