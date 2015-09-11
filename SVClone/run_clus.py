'''
Facilitates clustering of SVs
'''
import subprocess
import ipdb
import os
import pandas as pd
import scipy as sp
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import colorsys
import IPython
from collections import OrderedDict
from IPython.core.pylabtools import figsize

from . import cluster
from . import parameters as param
from . import cmd
from . import write_output

import numpy as np
from numpy import loadtxt
from numpy import zeros

def gen_new_colours(N):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def plot_cluster_hist(clusters,assignments,df,pl,pi,rlen,clus_out_dir,are_snvs=False):
    fig, axes = plt.subplots(1, 1, sharex=False, sharey=False, figsize=(12.5,4))
    RGB_tuples = gen_new_colours(len(clusters))
    
    sup,dep,cn_r,cn_v,mu_v,sides = [],[],[],[],[],[]
    Nvar = 0
    if are_snvs:
        clus_out_dir = '%s/snvs'%clus_out_dir
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)
        sup,ref,cn_r,cn_v,mu_v = cluster.get_snv_vals(df)
        dep = sup + ref
        av_cov = np.mean(dep)
        Nvar = len(sup)
        sides = np.zeros(Nvar,dtype=int)
    else:
        sup,dep,cn_r,cn_v,mu_v,sides,Nvar = cluster.get_sv_vals(df,rlen)

    for idx,clus in enumerate(clusters):
        clus_idx = np.array(assignments)==clus
        sup_clus = sup[clus_idx]
        dep_clus = dep[clus_idx]
        #axes[0].set_title("Clusters post-cluster merge: Cell fractions (raw VAFs purity-ploidy-adjusted)")
        #axes[0].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        axes.set_title("Raw VAFs")
        axes.hist(sup_clus/dep_clus,bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
    
    hist_plot = '%s/merged_cluster_hist'%clus_out_dir
    plt.savefig(hist_plot)

def plot_clusters(center_trace,clusters,assignments,snv_df,sv_df,pl,pi,rlen,clus_out_dir):
    fig, axes = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(12.5,8))

    RGB_tuples = gen_new_colours(len(clusters))

    axes[0].set_ylim([0,2])
    axes[0].set_title("Trace of $\phi_k$")

    for idx,clus in enumerate(clusters):
        axes[0].plot(center_trace[:, clus], label="trace of center %d" % clus, c=RGB_tuples[idx], lw=1)
        # print('phi_%d ~ %f' % (idx,np.mean(center_trace[2000:,clus])))

    leg = axes[0].legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)
    axes[1].set_title("Raw VAFs")
    #axes[2].set_title("Unmerged clusters: Cell fractions (raw VAFs purity-ploidy-adjusted)")
    
    sup,dep,cn_r,cn_v,mu_v,sides = [],[],[],[],[],[]
    Nvar = 0
    if are_snvs:
        clus_out_dir = '%s/snvs'%clus_out_dir
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)
        sup,ref,cn_r,cn_v,mu_v = cluster.get_snv_vals(df)
        dep = sup + ref
        av_cov = np.mean(dep)
        Nvar = len(sup)
        sides = np.zeros(Nvar,dtype=int)
    else:
        #sup,dep,cn_r,cn_v,mu_v,sides,Nvar = cluster.get_sv_vals(df,rlen)
        sup,dep,combos,sides,Nvar = cluster.get_sv_vals(df,rlen)

    for idx,clus in enumerate(clusters):
        clus_idx = np.array(assignments)==clus
        sup_clus = sup[clus_idx]
        dep_clus = dep[clus_idx]
        axes[1].hist(sup_clus/dep_clus,bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        #axes[2].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])

    plt.savefig('%s/cluster_trace_hist'%clus_out_dir)

def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def mean_confidence_interval(phi_trace, alpha=0.1):
    n = len(phi_trace)
    m, se = np.mean(phi_trace), np.std(phi_trace,ddof=1)/(np.sqrt(n))
    h = se * sp.stats.t._ppf((1+(1-alpha))/2., n-1)
    return round(m,4), round(m-h,4), round(m+h,4)

def merge_clusters(clus_out_dir,snv_df,sv_df,clus_info,clus_merged,clus_members,merged_ids,cparams,num_iters,burn,thin):
    '''
    cparams = [pi,rlen,insert,ploidy]
    '''
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
            mcmc = cluster.cluster(df.loc[new_members], cparams[0], cparams[1], \
                    cparams[2], cparams[3], num_iters, burn, thin, None, are_snvs, 1)
            trace = mcmc.trace("phi_k")[:]

            phis = mean_confidence_interval(trace)
            clus_merged.loc[idx] = np.array([ci.clus_id,new_size,phis[0],phis[1],phis[2]])            
            clus_members[idx] = new_members

            to_del.append(idx+1)
            merged_ids.append([int(ci.clus_id),int(cn.clus_id)])
            
            df_trace = pd.DataFrame(trace[:])
            df_trace.to_csv('%s/reclus%d_phi_trace.txt'%(clus_out_dir,int(ci.clus_id)),sep='\t',index=False)
            
            # plot new trace
            #assignments = [int(ci.clus_id)]*len(new_members)
            #tmp_trace = np.zeros((len(trace),int(ci.clus_id)+1))
            #tmp_trace[:,int(ci.clus_id)] = trace[:,0]
            #plot_clusters(tmp_trace,[int(ci.clus_id)],assignments,df.loc[new_members],cparams[3],cparams[0])

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
        return (clus_merged,clus_members,merged_ids)
    else: 
        new_df = pd.DataFrame(columns=clus_merged.columns,index=clus_merged.index)
        return merge_clusters(clus_out_dir,df,clus_merged,new_df,clus_members,merged_ids,cparams,num_iters,burn,thin,are_snvs)

def merge_results(clus_merged, merged_ids, df_probs, ccert):    
    # merge probability table
    to_drop = []
    df_probs_new = pd.DataFrame(df_probs,copy=True)
    
    for mid in merged_ids:
        clus1_vals = df_probs_new['cluster%d'%mid[0]].values
        clus2_vals = df_probs_new['cluster%d'%mid[1]].values
        
        df_probs_new['cluster%d'%mid[0]] = clus1_vals + clus2_vals
        to_drop.append(mid[1]) 

    to_drop = set(to_drop)
    for td in to_drop:
        df_probs_new = df_probs_new.drop('cluster%d'%td,1)
    
    # merge assignments in certainties table
    ccert_new = pd.DataFrame(ccert,copy=True)
    for mid in merged_ids:
        clus_max = ccert_new.most_likely_assignment.values
        ccert_new.loc[np.where(clus_max==mid[1])[0],'most_likely_assignment'] = mid[0]

    #update phis
    for idx,ci in clus_merged.iterrows():
        ccert_new.loc[ccert_new.most_likely_assignment==ci.clus_id,'average_ccf'] = ci.phi
        ccert_new.loc[ccert_new.most_likely_assignment==ci.clus_id,'90_perc_CI_lo'] = ci.phi_low_conf
        ccert_new.loc[ccert_new.most_likely_assignment==ci.clus_id,'90_perc_CI_hi'] = ci.phi_hi_conf

    return df_probs_new,ccert_new


def post_process_clusters(mcmc,sv_df,snv_df,merge_clusts,clus_out_dir,sample,pi,rlen):
    #npoints = len(df.spanning.values) if not are_snvs else len(df['var'].values)
    # assign points to highest probabiliy cluster
    npoints = len(snv_df) + len(sv_df)

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
    
    phis = np.array([mean_confidence_interval(center_trace[:,cid]) for cid in clus_info.clus_id.values])
    clus_info['phi'] = phis[:,0]
    clus_info['phi_low_conf'] = phis[:,1]
    clus_info['phi_hi_conf'] = phis[:,2]
    clus_info = clus_info.sort('phi',ascending=False)
    
    clus_ids = clus_info.clus_id.values
    clus_members = np.array([np.where(np.array(clus_max_prob)==i)[0] for i in clus_ids])

    col_names = map(lambda x: 'cluster'+str(x),clus_ids)
    df_probs = pd.DataFrame(clus_counts,dtype=float)[clus_ids].fillna(0)
    df_probs = df_probs.apply(lambda x: x/sum(x),axis=1)
    df_probs.columns = col_names

    # cluster certainty
    clus_max_df = pd.DataFrame(clus_max_prob,columns=['most_likely_assignment'])
    phi_cols = ["average_ccf","90_perc_CI_lo","90_perc_CI_hi"]
    phi_matrix = pd.DataFrame(phis[:],index=clus_ids,columns=phi_cols).loc[clus_max_prob]
    phi_matrix.index = range(len(phi_matrix))
    ccert = clus_max_df.join(phi_matrix)
    clus_info.index = range(len(clus_info))
    print(clus_info)

    if len(snv_df)>0 and len(sv_df)==0:
        # snvs only trace output
        clus_out_dir = '%s/snvs'%clus_out_dir        
    trace_out = 'premerge' if merge_clusts else ''
    trace_out = '%s/%s'%(clus_out_dir,trace_out)
    write_output.dump_trace(clus_info,center_trace,trace_out+'phi_trace.txt')
    write_output.dump_trace(clus_info,z_trace,trace_out+'z_trace.txt')
    
    if len(clus_info)>1 and merge_clusts: 
        #TODO: fix merging
        print('Merge clusters feature has been temporarily removed')

    # merge clusters
#    if len(clus_info)>1 and merge_clusts:        
#        clus_merged = pd.DataFrame(columns=clus_info.columns,index=clus_info.index)
#        clus_merged, clus_members, merged_ids  = merge_clusters(clus_out_dir,snv_df,sv_df,clus_info,clus_merged,\
#                clus_members,[],[pi,rlen,insert,ploidy],num_iters,burn,thin)
#        
#        if len(clus_merged)!=len(clus_info):
#            clus_info = clus_merged
#            df_probs, ccert = merge_results(clus_merged, merged_ids, df_probs, ccert)
#            if cmd.plot:
#                plot_cluster_hist(clus_merged.clus_id.values,ccert.most_likely_assignment.values,df,ploidy,pi,rlen,clus_out_dir,are_snvs)
    
    snv_probs = pd.DataFrame()
    snv_ccert = pd.DataFrame()
    snv_members = np.empty(0)
    if len(snv_df)>0:
        snv_pos = ['chrom','pos']
        snv_probs = df_probs.loc[:len(snv_df)-1]
        snv_probs = snv_df[snv_pos].join(snv_probs)
        
        snv_ccert = ccert.loc[:len(snv_df)-1]
        snv_ccert = snv_df[snv_pos].join(snv_ccert)

        snv_max_probs = np.array(clus_max_prob)[:len(snv_df)]
        snv_members = np.array([np.where(np.array(snv_max_probs)==i)[0] for i in clus_ids])

    sv_probs = pd.DataFrame()
    sv_ccert = pd.DataFrame()
    sv_members = np.empty(0)
    if len(sv_df)>0:
        lb = len(snv_df) if len(snv_df)>0 else 0
        
        sv_pos = ['bp1_chr','bp1_pos','bp2_chr','bp2_pos']
        sv_probs = df_probs.loc[lb:lb+len(sv_df)-1]
        sv_probs.index = sv_df.index
        sv_probs = sv_df[sv_pos].join(sv_probs)
        
        sv_ccert = ccert.loc[lb:lb+len(sv_df)-1]
        sv_ccert.index = sv_df.index
        sv_ccert = sv_df[sv_pos].join(sv_ccert)
        
        sv_max_probs = np.array(clus_max_prob)[:len(sv_df)]
        sv_members = np.array([np.where(np.array(sv_max_probs)==i)[0] for i in clus_ids])
    
    if len(clus_info) < 1:
        print("Warning! Could not converge on any major SV clusters. Skipping.\n")
    else:
        if len(snv_df)>0:
            write_output.write_out_files_snv(snv_df,clus_info,snv_members,
                    snv_probs,snv_ccert,clus_out_dir,sample,pi)
        if len(sv_df)>0:
            write_output.write_out_files(sv_df,clus_info,sv_members,
                    sv_probs,sv_ccert,clus_out_dir,sample,pi,rlen)
    
    # cluster plotting
    if cmd.plot:
        plot_clusters(center_trace,clus_idx,clus_max_prob,df,ploidy,pi,rlen,clus_out_dir,are_snvs)
    

#return clus_info,center_trace,z_trace,clus_members,df_probs,ccert

def infer_subclones(sample,sv_df,pi,rlen,insert,ploidy,out,n_runs,num_iters,burn,thin,beta,snv_df,merge_clusts,use_map,cocluster):
    clus_info,center_trace,ztrace,clus_members = None,None,None,None
    for i in range(n_runs):
        print("Cluster run: %d"%i)

        clus_out_dir = '%s/run%d'%(out,i)
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)    
        
        if cocluster:
            mcmc = cluster.cluster(sv_df,snv_df,pi,rlen,insert,ploidy,num_iters,burn,thin,beta,use_map)
            post_process_clusters(mcmc,sv_df,snv_df,merge_clusts,clus_out_dir,sample,pi,rlen)
        else:            
            if len(snv_df)>0:
                if not os.path.exists('%s/snvs'%clus_out_dir):
                    os.makedirs('%s/snvs'%clus_out_dir)
                mcmc = cluster.cluster(pd.DataFrame(),snv_df,pi,rlen,insert,ploidy,num_iters,burn,thin,beta,use_map)
                post_process_clusters(mcmc,pd.DataFrame(),snv_df,merge_clusts,clus_out_dir,sample,pi,rlen)
            if len(sv_df)>0:
                mcmc = cluster.cluster(sv_df,pd.DataFrame(),pi,rlen,insert,ploidy,num_iters,burn,thin,beta,use_map)
                post_process_clusters(mcmc,sv_df,pd.DataFrame(),merge_clusts,clus_out_dir,sample,pi,rlen)

