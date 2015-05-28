'''
Takes an SV matrix and attempts to build a phylogeny
'''
import subprocess
import ipdb
import sys
import os
import pandas as pd
import scipy as sp
import scipy.stats

from . import cluster
from . import parameters as param

import numpy as np
from numpy import loadtxt
from numpy import zeros

def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def mean_confidence_interval(phi_trace, alpha=0.1):
    n = len(phi_trace)
    m, se = np.mean(phi_trace), np.std(phi_trace,ddof=1)/(np.sqrt(n))
    h = se * sp.stats.t._ppf((1+(1-alpha))/2., n-1)
    return m, m-h, m+h

def merge_clusters(df,clus_info,clus_merged,clus_members,cparams,out,num_iters,burn,thin):
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
            trace = cluster.recluster(df.loc[new_members], cparams[0], cparams[1], cparams[2], num_iters, burn, thin)
            phis = mean_confidence_interval(trace)
            clus_members[idx] = new_members
            to_del.append(idx+1)
            df_trace = pd.DataFrame(np.transpose(trace[:]))
            df_trace.to_csv('%s/reclus%d_phi_trace.txt'%(out,int(ci.clus_id)),sep='\t',index=False)
            clus_merged.loc[idx] = np.array([ci.clus_id,new_size,phis[0],phis[1],phis[2]])
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
        return merge_clusters(df,clus_merged,new_df,clus_members,cparams,out,num_iters,burn,thin)

def run_clust(df,pi,rlen,insert,num_iters,burn,thin):
    mcmc = cluster.cluster(df,pi,rlen,insert,num_iters,burn,thin)
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
    # center_trace = center_trace[len(center_trace)/4:]
    # param.plot_clusters(center_trace,npoints,clus_idx)
    
    #phis = np.array([np.mean(center_trace[:,cid]) for cid in clus_info.clus_id.values])
    phis = np.array([mean_confidence_interval(center_trace[:,cid]) for cid in clus_info.clus_id.values])
    clus_info['phi'] = phis[:,0]
    clus_info['phi_low_conf'] = phis[:,1]
    clus_info['phi_hi_conf'] = phis[:,2]
    clus_info = clus_info.sort('phi',ascending=False)
    
    # filter clones below thresholds
    clus_info = clus_info[clus_info['phi'].values>(param.subclone_threshold/2)]
    clus_info = clus_info[sum(clus_info['size'].values)*param.subclone_sv_prop<clus_info['size'].values]
    clus_info.index = range(len(clus_info))
    
    clus_members = np.array([np.where(np.array(clus_max_prob)==i)[0] for i in clus_info.clus_id.values])

    # probability of being in each cluster
    probs = [map(lambda x: x[0]/x[1],zip(x,[float(sum(x))]*len(x))) for x in clus_counts]
    df_probs = pd.DataFrame(probs).fillna(0)
    df_pos = ['bp1_chr','bp1_pos','bp2_chr','bp2_pos']
    df_probs = df[df_pos].join(df_probs)
    # df_probs.columns= clusters etc.
    ipdb.set_trace()

    print(clus_info)
    return clus_info,center_trace,z_trace,clus_members,df_probs

def dump_trace(clus_info,center_trace,outf):
    traces = np.array([center_trace[:,cid] for cid in clus_info.clus_id.values])
    df_traces = pd.DataFrame(np.transpose(traces),columns=clus_info.clus_id)
    df_traces.to_csv(outf,sep='\t',index=False)

def write_out_files(df,clus_members,df_probs,clus_out_dir,sample):
    
    clus_info.to_csv('%s/clusters.txt'%(clus_out_dir),sep='\t',index=False)
    with open('%s/number_of_clusters.txt'%out,'w') as outf:
        outf.write("sample\tclusters\n")
        outf.write('%s\t%d\n'%(sample,len(clus_info)[0]))
    
    cn_dtype =  [('chr1','S50'),
                 ('pos1',int),
                 ('total_copynumber1',int),
                 ('no_chrs_bearing_mutation1',int),
                 ('chr2','S50'),
                 ('pos2',int),
                 ('total_copynumber2',int),
                 ('no_chrs_bearing_mutation2',int)]

    cmem = np.hstack(clus_members)
    cn_vect = np.empty((0,len(cmem)),dtype=cn_dtype)
    clus_svs = df.loc[cmem].copy()
    
    for idx,sv in clus_svs.iterrows():

        maj_cn1,min_cn1 = sv['bp1_maj_cnv'],sv['bp1_min_cnv']
        maj_cn2,min_cn2 = sv['bp2_maj_cnv'],sv['bp2_min_cnv']

        maj_cn1 = int(maj_cn1) if maj_cn1!="" else np.nan
        min_cn1 = int(min_cn1) if min_cn1!="" else np.nan
        maj_cn2 = int(maj_cn2) if maj_cn2!="" else np.nan
        min_cn2 = int(min_cn2) if min_cn2!="" else np.nan
        tot_cn1 = maj_cn1+min_cn1 if not np.isnan(maj_cn1) and not np.isnan(min_cn1) else np.nan
        tot_cn2 = maj_cn1+min_cn2 if not np.isnan(maj_cn2) and not np.isnan(min_cn2) else np.nan

        bp1_chr = str(sv['bp1_chr'])
        bp1_pos = int(sv['bp1_pos'])
        bp2_chr = str(sv['bp2_chr'])
        bp2_pos = int(sv['bp2_pos'])

        cn_new_row = np.array([(bp1_chr,bp1_pos,tot_cn1,min_cn1,bp2_chr,bp2_pos,tot_cn2,min_cn2)],dtype=cn_dtype)
        cn_vect = np.append(cn_vect,cn_new_row)

    pd.DataFrame(cn_vect).to_csv('%s/%s_copynumber.txt'%(clus_out_dir,sample),sep='\t',index=False)
    df_probs.to_csv('%s/%s_assignment_probability_table.txt'%(clus_out_dir,sample),sep='\t',index=False)

def infer_subclones(sample,df,pi,rlen,insert,out,n_runs,num_iters,burn,thin):

    clus_info,center_trace,ztrace,clus_members = None,None,None,None

    for i in range(n_runs):
        print("Cluster run: %d"%i)
        clus_info,center_trace,z_trace,clus_members,df_probs = run_clust(df,pi,rlen,insert,num_iters,burn,thin)
        sv_loss = 1-(sum(clus_info['size'])/float(len(df)))

        clus_out_dir = '%s/run%d'%(out,i)
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)

        with open("%s/warnings.txt"%clus_out_dir,'w') as warn:
            warn_str = ""
            if len(clus_info) < 1:
                warn_str = "Warning! Could not converge on any major clusters. Skipping.\n"
                warn.write(warn_str)
                print('\n'+warn_str)
                continue
            if sv_loss > param.tolerate_svloss:
                warn_str = "Warning! Lost %f of SVs.\n" % sv_loss
            if (clus_info.phi.values[0]*2.)>(1+param.subclone_diff):
                warn_str = warn_str + "Warning! Largest VAF cluster exceeds 0.5.\n"
            warn.write(warn_str)
            print('\n'+warn_str)
        
        dump_trace(clus_info,center_trace,'%s/phi_trace.txt'%clus_out_dir)
        dump_trace(clus_info,z_trace,'%s/z_trace.txt'%clus_out_dir)
        clus_info.to_csv('%s/clusters_premerged.txt'%(clus_out_dir),sep='\t',index=False)

#        clus_init = pd.DataFrame(columns=clus_info.columns,index=clus_info.index)
#        clus_info,clus_new_members = \
#            merge_clusters(df,clus_info,clus_init,clus_members.copy(),[pi,rlen,insert],clus_out_dir,num_iters,burn,thin)
#        
#        print("\nFiltered & merged clusters")
#        clus_info.clus_id = map(int,clus_info.clus_id.values)
#        clus_info['size'] = map(int,clus_info['size'].values)
#        print(clus_info)
        
        write_out_files(df,clus_info,clus_members,df_probs,clus_out_dir,sample)

