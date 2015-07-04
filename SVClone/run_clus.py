'''
Facilitates clustering of SVs
'''
import subprocess
import ipdb
import os
import pandas as pd
import scipy as sp
import scipy.stats
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

def plot_cluster_hist(clusters,assignments,df,pl,pi):
    fig, axes = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(12.5,8))
    RGB_tuples = gen_new_colours(len(clusters))

    for idx,clus in enumerate(clusters):
        t1 = df[np.array(assignments)==clus]
        s1 = t1.support.values
        n1 = map(np.mean,zip(t1.norm1.values,t1.norm2.values))
        axes[0].set_title("Clusters post-cluster merge: Cell fractions (raw VAFs purity-ploidy-adjusted)")
        axes[0].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        axes[1].set_title("Raw VAFs")
        axes[1].hist(s1/(n1+s1),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])

def plot_clusters(center_trace,clusters,assignments,df,pl,pi):
    fig, axes = plt.subplots(3, 1, sharex=False, sharey=False, figsize=(12.5,11))

    RGB_tuples = gen_new_colours(len(clusters))

    axes[0].set_ylim([0,1])
    axes[0].set_title("Trace of $\phi_k$")

    for idx,clus in enumerate(clusters):
        axes[0].plot(center_trace[:, clus], label="trace of center %d" % clus, c=RGB_tuples[idx], lw=1)
        # print('phi_%d ~ %f' % (idx,np.mean(center_trace[2000:,clus])))

    leg = axes[0].legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)

    for idx,clus in enumerate(clusters):
        t1 = df[np.array(assignments)==clus]
        s1 = t1.support.values
        n1 = map(np.mean,zip(t1.norm1.values,t1.norm2.values))
        axes[1].set_title("Unmerged clusters: Cell fractions (raw VAFs purity-ploidy-adjusted)")
        axes[1].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        #print 'mean cluster VAF: %f' % (np.mean(s1/(n1+s1)/pi))
        axes[2].set_title("Raw VAFs")
        axes[2].hist(s1/(n1+s1),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])

def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def mean_confidence_interval(phi_trace, alpha=0.1):
    n = len(phi_trace)
    m, se = np.mean(phi_trace), np.std(phi_trace,ddof=1)/(np.sqrt(n))
    h = se * sp.stats.t._ppf((1+(1-alpha))/2., n-1)
    return round(m,4), round(m-h,4), round(m+h,4)

def merge_clusters(clus_out_dir,df,clus_info,clus_merged,clus_members,merged_ids,cparams,num_iters,burn,thin,are_snvs=False):
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


def run_clust(clus_out_dir,df,pi,rlen,insert,ploidy,num_iters,burn,thin,beta,are_snvs=False):
    mcmc = cluster.cluster(df,pi,rlen,insert,ploidy,num_iters,burn,thin,beta,are_snvs)
    npoints = len(df.spanning.values) if not are_snvs else len(df['var'].values)

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
    
    phis = np.array([mean_confidence_interval(center_trace[:,cid]) for cid in clus_info.clus_id.values])
    clus_info['phi'] = phis[:,0]
    clus_info['phi_low_conf'] = phis[:,1]
    clus_info['phi_hi_conf'] = phis[:,2]
    clus_info = clus_info.sort('phi',ascending=False)
    
    clus_ids = clus_info.clus_id.values
    clus_members = np.array([np.where(np.array(clus_max_prob)==i)[0] for i in clus_ids])

    # TODO: reimplement cluster filtering based on size/floor threshold
    # cluster filtering
    #clus_info = clus_info[clus_info['phi'].values>(param.subclone_threshold/2)]
    #clus_info = clus_info[sum(clus_info['size'].values)*param.subclone_sv_prop<clus_info['size'].values]
    #clus_counts = [np.array(c)[clus_ids] for c in clus_counts]
    #clus_max_prob = np.array(clus_max_prob)[clus_ids]
    
    # probability of being in each cluster
    #probs = [np.array(x)[clus_info.clus_id.values] for x in probs]
    #clus_counts = [np.array(c)[clus_info.clus_id.values] for c in clus_counts]
    #probs = [map(lambda x: round(x[0]/x[1],4),zip(x,[float(sum(x))]*len(x))) for x in clus_counts]
    #df_probs = pd.DataFrame(probs).fillna(0)
    
    col_names = map(lambda x: 'cluster'+str(x),clus_ids)
    df_probs = pd.DataFrame(clus_counts,dtype=float)[clus_ids].fillna(0)
    df_probs = df_probs.apply(lambda x: x/sum(x),axis=1)
    df_probs.columns = col_names
    df_pos = []
    if are_snvs:
        df_pos = ['chrom','pos']
    else:
        df_pos = ['bp1_chr','bp1_pos','bp2_chr','bp2_pos']
    df_probs = df[df_pos].join(df_probs)

    # cluster certainty
    clus_max_df = pd.DataFrame(clus_max_prob,columns=['most_likely_assignment'])
    phi_cols = ["average_ccf","90_perc_CI_lo","90_perc_CI_hi"]
    phi_matrix = pd.DataFrame(phis[:],index=clus_ids,columns=phi_cols).loc[clus_max_prob]
    phi_matrix.index = range(len(phi_matrix))
    ccert = df[df_pos].join(clus_max_df).join(phi_matrix)

    clus_info.index = range(len(clus_info))
    print(clus_info)

    # cluster plotting
    if cmd.plot:
        plot_clusters(center_trace,clus_idx,clus_max_prob,df,ploidy,pi)
    write_output.dump_trace(clus_info,center_trace,'%s/premerge_phi_trace.txt'%clus_out_dir)
    write_output.dump_trace(clus_info,z_trace,'%s/premerge_z_trace.txt'%clus_out_dir)
    
    # merge clusters
    if len(clus_info)>1:        
        clus_merged = pd.DataFrame(columns=clus_info.columns,index=clus_info.index)
        clus_merged, clus_members, merged_ids  = merge_clusters(clus_out_dir,df,clus_info,clus_merged,\
                clus_members,[],[pi,rlen,insert,ploidy],num_iters,burn,thin,are_snvs)
        
        if len(clus_merged)!=len(clus_info):
            clus_info = clus_merged
            df_probs, ccert = merge_results(clus_merged, merged_ids, df_probs, ccert)
            if cmd.plot:
                plot_cluster_hist(clus_merged.clus_id.values,ccert.most_likely_assignment.values,df,ploidy,pi)
    
    return clus_info,center_trace,z_trace,clus_members,df_probs,ccert

def infer_subclones(sample,df,pi,rlen,insert,ploidy,out,n_runs,num_iters,burn,thin,beta,snv_df):
    if len(df) < 5:
        print("Less than 5 post-filtered SVs. Clustering not recommended for this sample. Exiting.")
        return None

    clus_info,center_trace,ztrace,clus_members = None,None,None,None
    for i in range(n_runs):
        print("Cluster run: %d"%i)

        clus_out_dir = '%s/run%d'%(out,i)
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)
    
        if len(snv_df)>0:
            print("Clustering SNVs...")
            clus_info,center_trace,z_trace,clus_members,df_probs,clus_cert = \
                run_clust(clus_out_dir,snv_df,pi,rlen,insert,ploidy,num_iters,burn,thin,beta,are_snvs=True)
            write_output.write_out_files_snv(snv_df,clus_info,clus_members,df_probs,clus_cert,clus_out_dir,sample,pi)

        print("Clustering SVs...")
        clus_info,center_trace,z_trace,clus_members,df_probs,clus_cert = \
                run_clust(clus_out_dir,df,pi,rlen,insert,ploidy,num_iters,burn,thin,beta)
        
        sv_loss = 1-(sum(clus_info['size'])/float(len(df)))

        if len(clus_info) < 1:
            print("Warning! Could not converge on any major clusters. Skipping.\n")
            continue
            #if sv_loss > param.tolerate_svloss:
            #    warn_str = "Warning! Lost %f of SVs.\n" % sv_loss
            #if (clus_info.phi.values[0]*2.)>(1+param.subclone_diff):
            #    warn_str = warn_str + "Warning! Largest VAF cluster exceeds 0.5.\n"
        
        write_output.write_out_files(df,clus_info,clus_members,df_probs,clus_cert,clus_out_dir,sample,pi,rlen)

