'''
Facilitates clustering of SVs
'''
from __future__ import print_function

import subprocess
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import colorsys
import IPython
import multiprocessing
import time
import shutil
import pymc as pm
import random
import ConfigParser

from distutils.dir_util import copy_tree
from collections import OrderedDict
from IPython.core.pylabtools import figsize
from pymc.utils import hpd

from . import cluster
from . import load_data
from . import write_output
from SVprocess import svp_load_data as svp_load

import numpy as np
from numpy import loadtxt
from numpy import zeros

def gen_new_colours(N):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def plot_clusters(trace, clusters, assignments, sup, dep, clus_out_dir, cparams):

    center_trace = trace("phi_k")[:]
    z_trace = trace('z')[:]
    phi_limit = cparams['phi_limit']
    maxclus = max([max(z) for z in z_trace])

    alpha_trace = []
    try:
        alpha_trace = trace('alpha')[:]
    except KeyError:
        pass

    fig, axes = plt.subplots(4, 1, sharex=False, sharey=False, figsize=(12.5,12))
    if len(alpha_trace) > 0:
        fig, axes = plt.subplots(5, 1, sharex=False, sharey=False, figsize=(12.5,15))
        axes[4].set_title("Alpha trace")
        axes[4].plot(range(len(alpha_trace)), alpha_trace, lw=1)

    RGB_tuples = gen_new_colours(len(clusters))

    axes[0].set_ylim([0, phi_limit + 0.1])
    axes[0].set_xlim([0, len(center_trace)])
    axes[0].set_title("Trace of $\phi_k$")

    axes[1].set_ylim([0, phi_limit + 0.1])
    axes[1].set_xlim([0, len(center_trace)])
    axes[1].set_title("Adjusted trace of $\phi_k$")

    axes[2].set_ylim([-1, maxclus+1])
    axes[2].set_xlim([0, len(z_trace)])
    axes[2].set_title("Trace of z categorical")

    axes[3].set_title("Raw VAFs")

    center_trace_adj = get_adjusted_phi_trace(center_trace, clusters)
    for idx,clus in enumerate(clusters):
        axes[0].plot(center_trace[:, clus], label="trace of center %d" % clus, c=RGB_tuples[idx], lw=1)
        axes[1].plot(center_trace_adj[:, clus], label="adjusted trace of center %d" % clus, c=RGB_tuples[idx], lw=1)

    leg = axes[0].legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)

    for i in range(len(z_trace[0])):
        idx = int(np.where(assignments[i]==clusters)[0])
        axes[2].plot(z_trace[:, i], label="", c=RGB_tuples[idx], lw=1, alpha=0.2)

    for idx,clus in enumerate(clusters):
        clus_idx = np.array(assignments)==clus
        sup_clus = sup[clus_idx]
        dep_clus = dep[clus_idx]
        axes[3].hist(sup_clus/dep_clus,bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])

    fig.tight_layout()
    plt.savefig('%s/cluster_trace_hist'%clus_out_dir)

def index_max(values):
    return max(xrange(len(values)),key=values.__getitem__)

def mean_confidence_interval(phi_trace, alpha):
    n = len(phi_trace)
    m = np.mean(phi_trace)
    #m, se = np.mean(phi_trace), np.std(phi_trace,ddof=1)/(np.sqrt(n))
    #h = se * sp.stats.t._ppf((1+(1-alpha))/2., n-1)
    phi_hpd = hpd(phi_trace, alpha)
    return round(m,4), round(phi_hpd[0],4), round(phi_hpd[1],4) #round(m-h,4), round(m+h,4)

def merge_clusters(clus_out_dir,clus_info,clus_merged,clus_members,merged_ids,sup,dep,norm,cn_states,sparams,cparams):
    clus_limit = cparams['clus_limit']
    phi_limit = cparams['phi_limit']
    subclone_diff = cparams['subclone_diff']
    norm = np.array(norm)

    if len(clus_info)==1:
        return (clus_info,clus_members)

    clus_info = clus_info.copy()
    clus_info = clus_info.sort_values('phi', ascending=False)
    clus_info.index = range(0,len(clus_info))
    to_del = []
    for idx,ci in clus_info.iterrows():
        if idx+1 >= len(clus_info):
            clus_merged.loc[idx] = clus_info.loc[idx]
            break
        cn = clus_info.loc[idx+1]
        # low CI boundary of current cluster, high CI boundary of next cluster
        ci_lo, cn_hi = ci[3], cn[4]
        if cn_hi >= (ci_lo - subclone_diff):
            print("\nReclustering similar clusters...")
            new_size = ci['size'] + cn['size']

            new_members = np.concatenate([clus_members[idx],clus_members[idx+1]])
            mcmc, map_ = cluster.cluster(sup[new_members], dep[new_members],
                                         cn_states[new_members], len(new_members),
                                         sparams,cparams, phi_limit, norm[new_members], recluster=True)
            trace = mcmc.trace("phi_k")[:]

            phis = mean_confidence_interval(trace,cparams['hpd_alpha'])
            clus_merged.loc[idx] = np.array([ci.clus_id,new_size,phis[0],phis[1],phis[2],phis[0],phis[1],phis[2]])
            clus_members[idx] = new_members

            to_del.append(idx+1)
            merged_ids.append([int(ci.clus_id),int(cn.clus_id)])

            df_trace = pd.DataFrame(trace[:])
            df_trace.to_csv('%s/reclus%d_phi_trace.txt'%(clus_out_dir,int(ci.clus_id)),sep='\t',index=False)

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
        clus_merged = clus_merged.sort_values('phi',ascending=False)
        clus_merged.index = range(len(clus_merged))

        print("\nMerged clusters")
        print('Input')
        print(clus_info)
        print('Merged')
        print(clus_merged)

    if len(clus_merged) == 1 or (len(clus_info) == len(clus_merged)):
        return (clus_merged,clus_members,merged_ids)
    else:
        new_df = pd.DataFrame(columns=clus_merged.columns,index=clus_merged.index)
        return merge_clusters(clus_out_dir,clus_merged,new_df,clus_members,merged_ids,sup,dep,norm,cn_states,sparams,cparams)

def merge_results(clus_merged, merged_ids, df_probs, ccert):
    # merge probability table
    to_drop = []
    df_probs_new = pd.DataFrame(df_probs,copy=True)
    hpd_lo = ccert.columns.values[2]
    hpd_hi = ccert.columns.values[3]

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
        ccert_new.loc[ccert_new.most_likely_assignment==ci.clus_id, 'average_ccf'] = ci['phi']
        ccert_new.loc[ccert_new.most_likely_assignment==ci.clus_id, hpd_lo] = ci[hpd_lo]
        ccert_new.loc[ccert_new.most_likely_assignment==ci.clus_id, hpd_hi] = ci[hpd_hi]

    return df_probs_new,ccert_new

def get_adjusted_phi_trace(center_trace, clus_idx):
    center_trace_adj = center_trace.copy()
    if len(clus_idx) > 1:
        for i in range(len(center_trace)):
            ph = center_trace[i][clus_idx]
            ranks = ph.argsort()[::-1]
            for idx,clus in enumerate(clus_idx):
                center_trace_adj[i][clus] = ph[ranks[idx]]
        return(center_trace_adj)
    else:
        return(center_trace_adj)

def get_adjusted_phis(clus_info, center_trace, cparams):
    '''
    Fixes potential label-switching problems by re-ordering phi traces from
    smallest to largest then assigning phis in the order of the unadjusted phis
    '''
    hpd_alpha     = cparams['hpd_alpha']

    clus_idx = clus_info.clus_id.values
    center_trace_adj = get_adjusted_phi_trace(center_trace, clus_idx)
    phis = np.array([mean_confidence_interval(center_trace[:,cid],hpd_alpha) for cid in clus_idx])

    if len(clus_idx) > 1:
        phis_adj = np.array([mean_confidence_interval(center_trace_adj[:,cid],hpd_alpha) for cid in clus_idx])
        phis_sort = np.argsort(phis[:,0])[::-1]

        for i in range(len(phis_sort)):
            on_idx = phis_sort[i]
            phis[on_idx] = phis_adj[i]

        return(phis)
    else:
        return(phis)

def get_per_variant_phi(z_trace, phi_trace):

    z_phi = np.array(z_trace.copy(), dtype=float)
    for i in range(len(z_trace)):
        for j in range(len(z_trace[0])):
            clus = int(z_trace[i,j])
            z_phi[i,j] = phi_trace[i, clus]

    return z_phi.mean(axis=0)

def post_process_clusters(mcmc,sv_df,snv_df,clus_out_dir,sup,dep,norm,cn_states,sparams,cparams,output_params,map_):

    merge_clusts  = cparams['merge_clusts']
    subclone_diff = cparams['subclone_diff']
    phi_limit     = cparams['phi_limit']
    merge_clusts  = cparams['merge_clusts']
    cnv_pval      = cparams['clonal_cnv_pval']
    hpd_alpha     = cparams['hpd_alpha']
    adjust_phis   = cparams['adjust_phis']
    clus_penalty  = output_params['cluster_penalty']
    smc_het       = output_params['smc_het']
    plot          = output_params['plot']

    try:
        sv_df = sv_df[sv_df.classification.values!='SIMU_SV']
    except AttributeError:
        pass
    npoints = len(snv_df) + len(sv_df)
    sup, dep, norm, cn_states = sup[:npoints], dep[:npoints], norm[:npoints], cn_states[:npoints]

    z_trace = mcmc.trace('z')[:]

    # assign points to highest probability cluster
    clus_counts = [np.bincount(z_trace[:,i]) for i in range(npoints)]
    clus_max_prob = [index_max(c) for c in clus_counts]
    clus_mp_counts = np.bincount(clus_max_prob)
    clus_idx = np.nonzero(clus_mp_counts)[0]
    clus_mp_counts = clus_mp_counts[clus_idx]

    # cluster distribution
    clus_info = pd.DataFrame(clus_idx, columns=['clus_id'])
    clus_info['size'] = clus_mp_counts

    if len(clus_info) < 1:
        print("Warning! Could not converge on any major SV clusters. Skipping.\n")
        return None

    center_trace = mcmc.trace("phi_k")[:]
    phis = np.array([mean_confidence_interval(center_trace[:,cid], hpd_alpha) for cid in clus_idx])
    original_phis = phis.copy()
    adjusted_phis = get_adjusted_phis(clus_info, center_trace, cparams)

    hpd_lo = '_'.join([str(int(100-(100*hpd_alpha))), 'HPD', 'lo'])
    hpd_hi = '_'.join([str(int(100-(100*hpd_alpha))), 'HPD', 'hi'])

    phis = adjusted_phis if adjust_phis else phis
    clus_info['phi']  = phis[:,0]
    clus_info[hpd_lo] = phis[:,1]
    clus_info[hpd_hi] = phis[:,2]

    if adjust_phis:
        clus_info['phi_unadjusted']  = original_phis[:,0]
        clus_info['%s_unadjusted'%hpd_lo] = original_phis[:,1]
        clus_info['%s_unadjusted'%hpd_hi] = original_phis[:,2]
    else:
        clus_info['phi_adjusted']  = adjusted_phis[:,0]
        clus_info['%s_adjusted'%hpd_lo] = adjusted_phis[:,1]
        clus_info['%s_adjusted'%hpd_hi] = adjusted_phis[:,2]

    clus_ids = clus_info.clus_id.values
    clus_members = np.array([np.where(np.array(clus_max_prob)==i)[0] for i in clus_ids])

    col_names = map(lambda x: 'cluster'+str(x),clus_ids)
    df_probs = pd.DataFrame(clus_counts,dtype=float)[clus_ids].fillna(0)
    df_probs = df_probs.apply(lambda x: x/sum(x),axis=1)
    df_probs.columns = col_names

    # cluster certainty
    clus_max_df = pd.DataFrame(clus_max_prob,columns=['most_likely_assignment'])
    phi_cols = ["average_ccf", hpd_lo, hpd_hi]
    phi_matrix = pd.DataFrame(phis[:],index=clus_ids,columns=phi_cols).loc[clus_max_prob]
    phi_matrix.index = range(len(phi_matrix))
    ccert = clus_max_df.join(phi_matrix)
    clus_info.index = range(len(clus_info))

    print('\n\n')
    print(clus_info[['clus_id', 'size', 'phi']])
    print('Compiling and writing output...')

    dump_out_dir = clus_out_dir
    if len(snv_df)>0 and len(sv_df)==0:
        # snvs only trace output
        dump_out_dir = '%s/snvs'%clus_out_dir
    trace_out = '%s/' % (dump_out_dir)
    write_output.dump_trace(center_trace, trace_out+'phi_trace.txt')
    write_output.dump_trace(z_trace, trace_out+'z_trace.txt')

    try:
        alpha_trace = mcmc.trace('alpha')[:]
        write_output.dump_trace(alpha_trace, trace_out+'alpha_trace.txt')
    except KeyError:
        pass

    # cluster plotting
    if plot:
        plot_clusters(mcmc.trace, clus_idx, clus_max_prob, sup, dep, clus_out_dir, cparams)

    # merge clusters
    if len(clus_info)>1 and merge_clusts:
        clus_merged = pd.DataFrame(columns=clus_info.columns,index=clus_info.index)
        clus_merged, clus_members, merged_ids  = merge_clusters(clus_out_dir,clus_info,clus_merged,\
                clus_members,[],sup,dep,norm,cn_states,sparams,cparams)

        if len(clus_merged)!=len(clus_info):
            clus_info = clus_merged
            df_probs, ccert = merge_results(clus_merged, merged_ids, df_probs, ccert)

    snv_probs = pd.DataFrame()
    snv_ccert = pd.DataFrame()
    snv_members = np.empty(0)

    z_phi = get_per_variant_phi(z_trace, center_trace)

    # compile run fit statistics
    run_fit = pd.DataFrame()
    if map_ is not None:
        nclus = len(clus_info)
        # bic = -2 * map_.lnL + (1 + npoints + nclus * 2) + (nclus * clus_penalty) * np.log(npoints)
        phis = ccert.average_ccf.values
        cns, pvs = cluster.get_most_likely_cn_states(cn_states, sup, dep, phis, sparams['pi'], cnv_pval, norm)
        lls = []
        for si, di, pvi in zip(sup, dep, pvs):
            lls.append(pm.binomial_like(si, di, pvi))
        svc_ic = -2 * np.sum(lls) + (npoints + nclus * clus_penalty) * np.log(npoints)

        run_fit = pd.DataFrame([['svc_IC', svc_ic], ['BIC', map_.BIC], ['AIC', map_.AIC], ['AICc', map_.AICc],
                                ['lnL', map_.lnL], ['logp', map_.logp], ['logp_at_max', map_.logp_at_max],
                                ['param_len', map_.len], ['data_len', map_.data_len]])

    if len(snv_df)>0:
        snv_pos = ['chrom','pos']
        snv_probs = df_probs.loc[:len(snv_df)-1]
        snv_probs = snv_df[snv_pos].join(snv_probs)

        snv_ccert = ccert.loc[:len(snv_df)-1]
        snv_ccert = snv_df[snv_pos].join(snv_ccert)

        snv_max_probs = np.array(clus_max_prob)[:len(snv_df)]
        snv_members = np.array([np.where(np.array(snv_max_probs)==i)[0] for i in clus_ids])

        snv_sup       = sup[:len(snv_df)]
        snv_dep       = dep[:len(snv_df)]
        snv_norm      = norm[:len(snv_df)]
        snv_cn_states = cn_states[:len(snv_df)]
        snv_z_phi     = z_phi[:len(snv_df)]
        write_output.write_out_files(snv_df,clus_info.copy(),snv_members,
                snv_probs,snv_ccert,clus_out_dir,sparams['sample'],sparams['pi'],snv_sup,
                snv_dep,snv_norm,snv_cn_states,run_fit,smc_het,cnv_pval,snv_z_phi,are_snvs=True)

    sv_probs = pd.DataFrame()
    sv_ccert = pd.DataFrame()
    sv_members = np.empty(0)
    if len(sv_df)>0:
        lb = len(snv_df) if len(snv_df)>0 else 0

        sv_pos = ['chr1','pos1','dir1','chr2','pos2','dir2']
        sv_probs = df_probs.loc[lb:lb+len(sv_df)-1]
        sv_probs.index = sv_df.index
        sv_probs = sv_df[sv_pos].join(sv_probs)

        sv_ccert = ccert.loc[lb:lb+len(sv_df)-1]
        sv_ccert.index = sv_df.index
        sv_ccert = sv_df[sv_pos].join(sv_ccert)

        sv_max_probs = np.array(clus_max_prob)[:len(sv_df)]
        sv_members = np.array([np.where(np.array(sv_max_probs)==i)[0] for i in clus_ids])

        sv_sup       = sup[lb:lb+len(sv_df)]
        sv_dep       = dep[lb:lb+len(sv_df)]
        sv_norm      = norm[lb:lb+len(sv_df)]
        sv_cn_states = cn_states[lb:lb+len(sv_df)]
        sv_z_phi     = z_phi[lb:lb+len(sv_df)]
        write_output.write_out_files(sv_df,clus_info.copy(),sv_members,
                    sv_probs,sv_ccert,clus_out_dir,sparams['sample'],sparams['pi'],sv_sup,
                    sv_dep,sv_norm,sv_cn_states,run_fit,smc_het,cnv_pval,sv_z_phi)

def cluster_and_process(sv_df, snv_df, run, out_dir, sample_params, cluster_params, output_params, seeds):
    np.random.seed(seeds[run])
    print('Random seed for run %d is %d' % (run, seeds[run]))

    clus_out_dir = '%s/run%d'%(out_dir, run)
    if not os.path.exists(clus_out_dir):
        os.makedirs(clus_out_dir)

    Nvar, sup, dep, cn_states, norm = None, None, None, None, None
    if cluster_params['cocluster'] and len(sv_df)>0 and len(snv_df)>0:
        # coclustering
        sup, dep, cn_states, Nvar, norm = load_data.get_snv_vals(snv_df, cluster_params)
        sv_sup, sv_dep, sv_cn_states, sv_Nvar, sv_norm = load_data.get_sv_vals(sv_df,
                                                cluster_params['adjusted'], cluster_params)
        sup = np.append(sup, sv_sup)
        dep = np.append(dep, sv_dep)
        norm = np.append(norm, sv_norm)
        cn_states = pd.concat([pd.DataFrame(cn_states),pd.DataFrame(sv_cn_states)])[0].values
        Nvar = Nvar + sv_Nvar
        mcmc, map_ = cluster.cluster(sup, dep, cn_states, Nvar, sample_params,
                                     cluster_params, cluster_params['phi_limit'], norm)
        post_process_clusters(mcmc, sv_df, snv_df, clus_out_dir, sup, dep, norm, cn_states,
                          sample_params, cluster_params, output_params, map_)

    elif len(sv_df) > 0 or len(snv_df) > 0:
        # no coclustering
        if len(snv_df) > 0:
            if not os.path.exists('%s/snvs'%clus_out_dir):
                os.makedirs('%s/snvs'%clus_out_dir)
            sup, dep, cn_states, Nvar, norm = load_data.get_snv_vals(snv_df, cluster_params)
            mcmc, map_ = cluster.cluster(sup, dep, cn_states, Nvar, sample_params,
                                         cluster_params, cluster_params['phi_limit'], norm)
            post_process_clusters(mcmc, pd.DataFrame(), snv_df, clus_out_dir, sup, dep, norm,
                                  cn_states,sample_params, cluster_params, output_params, map_)
        if len(sv_df) > 0:
            sup, dep, cn_states, Nvar, norm = load_data.get_sv_vals(sv_df,
                                                cluster_params['adjusted'], cluster_params)
            mcmc, map_ = cluster.cluster(sup, dep, cn_states, Nvar, sample_params,
                                         cluster_params, cluster_params['phi_limit'], norm)
            post_process_clusters(mcmc, sv_df, pd.DataFrame(), clus_out_dir, sup, dep, norm,
                                  cn_states, sample_params, cluster_params, output_params, map_)

    else:
        raise ValueError("No valid variants to cluster!")

def pick_best_run(n_runs, out, sample, ccf_reject, cocluster, fit_metric, cluster_penalty, are_snvs=False):
    snv_dir = 'snvs/' if are_snvs else ''

    ics = []
    for run in range(n_runs):
        fit_file = '%s/run%d/%s%s_fit.txt' % (out, run, snv_dir, sample)
        fit = pd.read_csv(fit_file, delimiter='\t', dtype=None, header=None)
        try:
            ic_val = fit.loc[fit[0]==fit_metric].values[0][1]
            ics.append(ic_val)
        except IndexError:
            raise ValueError('Invalid fit metric specified or fit metric does not exist in run.')
    ics = np.array(ics)

    min_ic = -1
    ic_sort = np.argsort(ics)
    for idx in range(n_runs):
        min_ic = ic_sort[idx]
        struct_file = '%s/run%d/%s%s_subclonal_structure.txt' % (out, min_ic, snv_dir, sample)
        clus_struct = pd.read_csv(struct_file, delimiter='\t', dtype=None, header=0)

        if len(clus_struct) > 1:
            var_props = clus_struct.n_variants.map(float).values/sum(clus_struct.n_variants)
            clus_struct = clus_struct[var_props > 0.01] #filter out very small clusters
            if len(clus_struct) > 1:
                break

        if clus_struct.CCF[0] > ccf_reject:
            print('Successfully finished clustering! Best run:')
            print(clus_struct[['cluster', 'n_variants', 'proportion', 'CCF']])
            break
        else:
            min_ic = -1

    if min_ic == -1:
        print('No optimal run found! Consider more runs and/or more iterations.')
    else:
        if cocluster:
            print('Selecting run %d as best run for SVs & SNVs' % min_ic)
            best_run = '%s/run%d' % (out, min_ic)
            best_run_dest = '%s/best_run_coclus' % out
            copy_tree(best_run, best_run_dest)

            f = open('%s/run%s.txt' % (best_run_dest, min_ic), 'w')
            f.close()

        elif are_snvs:
            print('Selecting run %d as best run for SNVs' % min_ic)
            best_run = '%s/run%d/snvs' % (out, min_ic)
            best_run_dest = '%s/best_run_snvs' % out
            copy_tree(best_run, best_run_dest)

            shutil.copyfile('%s/run%d/snvs/phi_trace.txt.gz' % (out, min_ic), '%s/phi_trace.txt.gz' % best_run_dest)
            shutil.copyfile('%s/run%d/snvs/z_trace.txt.gz' % (out, min_ic), '%s/z_trace.txt.gz' % best_run_dest)
            shutil.copyfile('%s/run%d/cluster_trace_hist.png' % (out, min_ic),
                            '%s/cluster_trace_hist.png' % best_run_dest)

            f = open('%s/run%s.txt' % (best_run_dest, min_ic), 'w')
            f.close()
        else:
            print('Selecting run %d as best run for SVs' % min_ic)
            best_run = '%s/run%d' % (out, min_ic)
            best_run_dest = '%s/best_run_svs' % out
            copy_tree(best_run, best_run_dest)

            snv_folder = '%s/best_run_svs/snvs' % out
            if os.path.exists(snv_folder):
                shutil.rmtree(snv_folder)

            f = open('%s/run%s.txt' % (best_run_dest, min_ic), 'w')
            f.close()

def simu_sv(sv):
    simu_sv = sv.copy()
    simu_sv.classification = 'SIMU_SV'
    simu_sv[8:30] = 0
    simu_sv.raw_mean_vaf = 0.0
    simu_sv.adjusted_support = np.random.binomial(sv.adjusted_depth, sv.adjusted_vaf)
    simu_sv.adjusted_vaf = float(simu_sv.adjusted_support) / simu_sv.adjusted_depth
    return simu_sv

def get_seeds(seeds, n_runs):
    '''
    splits seed argument list or if no list,
    generate a new list of random seeds
    '''
    try:
        if seeds =='':
            seeds = [np.random.randint(0,10000) for i in range(n_runs)]
        else:
            seeds = [int(s) for s in seeds.split(',')]
    except ValueError:
        seeds = [np.random.randint(0,10000) for i in range(n_runs)]

    return(seeds)

def subsample_snvs(snv_df, subsample, run, ss_seeds, sample, out):
    print('Subsampling %d SNVs for run%d' % (subsample, run))
    print('Random seed for sampling is %d' % ss_seeds[run])

    np.random.seed(ss_seeds[run])
    keep = np.random.choice(snv_df.index, size=subsample, replace=False)
    snv_run_df = snv_df.loc[keep]

    snv_run_df.index = range(len(snv_run_df)) #re-index
    snv_outname = '%s/%s_filtered_snvs_%d_subsampled_run%d.tsv' % (out, sample, subsample, run)
    snv_run_df.to_csv(snv_outname, sep='\t', index=False, na_rep='')

    return(snv_run_df)

def select_copynumber(cn_state):
    '''
    Return the major, minor and subclonal fraction
    for clonal or subclonal copy-number states
    '''
    try:
        cn_state = [cn.split(',') for cn in cn_state.split('|')]
        cn_state = [float(cn) for cn_side in cn_state for cn in cn_side]

        if len(cn_state) < 4:
            # extend CN array with zeroes if clonal
            cn_state.extend([0.0] * 3)
            return cn_state
        else:
            return cn_state
    except:
        return [0] * 6

def format_snvs_for_ccube(df, sparams, cparams, cc_file):
    '''
    prepare ccube input for SNVs
    '''
    sup, dep, cn_states, Nvar, norm_cn = load_data.get_snv_vals(df, cparams)
    mutation_id = ['%s_%d' % (c, p) for c, p in zip(df.chrom, df.pos)]

    cn_states = df.gtype.map(lambda x: select_copynumber(x))
    cn_states = pd.DataFrame.from_records(list(cn_states.values))
    cc_df = pd.DataFrame({'mutation_id': mutation_id,
                         'purity': sparams['pi'],
                         'normal_cn': norm_cn,
                         'var_counts': sup,
                         'ref_counts': dep - sup,
                         'total_counts': dep,
                         'vaf': sup / dep,
                         'subclonal_cn': cn_states[2] < 1,
                         'major_cn_sub1': cn_states[0],
                         'minor_cn_sub1': cn_states[1],
                         'total_cn_sub1': cn_states[0] + cn_states[1],
                         'frac_cn_sub1': cn_states[2],
                         'major_cn_sub2': cn_states[3],
                         'minor_cn_sub2': cn_states[4],
                         'total_cn_sub2': cn_states[3] + cn_states[4],
                         'frac_cn_sub2': cn_states[5]})

    cc_df.to_csv(cc_file, sep='\t', index=False)

def format_svs_for_ccube(df, sparams, cparams, cc_file, gain_loss_classes):
    '''
    prepare ccube input for SVs
    '''
    dl = len(df)
    df = df[np.logical_and(df.gtype1!='', df.gtype2!='')]
    if len(df) != dl:
        print('Filtered out %d SVs with missing copy-number states' % (dl - len(df)))

    adjust_factor = 1. - (float(sparams['pi']) / sparams['ploidy'])
    sup, dep, cn_states, Nvar, norm_cn = load_data.get_sv_vals(df, cparams)

    sv_classes = df.classification.values
    dups = np.array([sv_class in gain_loss_classes['dna_gain_class'] for idx, sv_class in enumerate(sv_classes)])
    norm1, norm2 = df.norm1, df.norm2
    norm1[dups] = norm1[dups] * adjust_factor
    norm2[dups] = norm2[dups] * adjust_factor

    mutation_id = ['%s:%d:%s_%s:%d:%s' % (c1, p1, d1, c2, p2, d2) for c1, p1, d1, c2, p2, d2 in \
                    zip(df.chr1.values, df.pos1.values, df.dir1.values, df.chr2.values, df.pos2.values, df.dir2.values)]

    cn_states1 = df.gtype1.map(lambda x: select_copynumber(x))
    cn_states2 = df.gtype2.map(lambda x: select_copynumber(x))
    cn_states1 = pd.DataFrame.from_records(list(cn_states1.values))
    cn_states2 = pd.DataFrame.from_records(list(cn_states2.values))

    cc_df = pd.DataFrame({'id': df.ID.values,
                         'mutation_id': mutation_id,
                         'purity': sparams['pi'],
                         'normal_cn': norm_cn,
                         'var_counts1': sup,
                         'ref_counts1': norm1.values,
                         'total_counts1': sup + norm1.values,
                         'vaf1': sup / (sup + norm1.values),
                         'subclonal_cn1': cn_states1[2] < 1,
                         'major_cn1_sub1': cn_states1[0],
                         'minor_cn1_sub1': cn_states1[1],
                         'total_cn1_sub1': cn_states1[0] + cn_states1[1],
                         'frac_cn1_sub1': cn_states1[2],
                         'major_cn1_sub2': cn_states1[3],
                         'minor_cn1_sub2': cn_states1[4],
                         'total_cn1_sub2': cn_states1[3] + cn_states1[4],
                         'frac_cn1_sub2': cn_states1[5],
                         'var_counts2': sup,
                         'ref_counts2': norm2.values,
                         'total_counts2': sup + norm2.values,
                         'vaf2': sup / (sup + norm2.values),
                         'subclonal_cn2': cn_states2[2] < 1,
                         'major_cn2_sub1': cn_states2[0],
                         'minor_cn2_sub1': cn_states2[1],
                         'total_cn2_sub1': cn_states2[0] + cn_states2[1],
                         'frac_cn2_sub1': cn_states2[2],
                         'major_cn2_sub2': cn_states2[3],
                         'minor_cn2_sub2': cn_states2[4],
                         'total_cn2_sub2': cn_states2[3] + cn_states2[4],
                         'frac_cn2_sub2': cn_states2[5]})

    cc_df.to_csv(cc_file, sep='\t', index=False)

def run_clustering(args):
    sample          = args.sample
    cfg             = args.cfg
    out             = sample if args.out == "" else args.out
    pp_file         = args.pp_file
    param_file      = args.param_file
    snv_file        = args.snv_file
    sv_file         = args.sv_file
    seeds           = args.seeds
    XX              = args.XX
    XY              = args.XY
    subsample       = args.subsample
    ss_seeds        = args.ss_seeds

    if out!='' and not os.path.exists(out):
        os.makedirs(out)

    sample_params, cluster_params, output_params = \
                   load_data.get_params_cluster_step(sample, cfg, out, pp_file, param_file, XX, XY)
    n_runs = cluster_params['n_runs']

    seeds = get_seeds(seeds, n_runs)
    ss_seeds = get_seeds(ss_seeds, n_runs)

    sv_df       = pd.DataFrame()
    snv_df      = pd.DataFrame()
    sv_file     = sv_file if sv_file != '' else '%s/%s_filtered_svs.tsv' % (out,sample)
    snv_file    = snv_file if snv_file != '' else '%s/%s_filtered_snvs.tsv' % (out,sample)

    if not os.path.exists(snv_file) and not os.path.exists(sv_file):
        raise OSError('Neither SNV nor SV input files were found!')

    if os.path.exists(sv_file):
        sv_df = pd.read_csv(sv_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        sv_df = pd.DataFrame(sv_df).fillna('')

    if os.path.exists(snv_file):
        snv_df = pd.read_csv(snv_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        snv_df = pd.DataFrame(snv_df).fillna('')

    if (len(sv_df) == 0 or len(snv_df) ==0) and cluster_params['cocluster']:
        cluster_params['cocluster'] = False

    clus_info, center_trace, ztrace, clus_members = None, None, None, None

    sv_to_sim = cluster_params['sv_to_sim']
    if sv_to_sim > 0:
        print('Simulating %d SV...' % (len(sv_df) * sv_to_sim))
        nvar = len(sv_df)
        for i in range(nvar):
            for j in range(sv_to_sim):
                sv_df.loc[len(sv_df)] = simu_sv(sv_df.loc[i])

    threads = cluster_params['threads']
    use_map = cluster_params['use_map']
    ccf_reject = cluster_params['ccf_reject']
    cocluster = cluster_params['cocluster']
    fit_metric = output_params['fit_metric']
    cluster_penalty = output_params['cluster_penalty']

    if cocluster:
        print("Coclustering %d SVs & %d SNVs..." % (len(sv_df), len(snv_df)))
    else:
        if len(snv_df) > 0:
            print("Clustering %d SNVs..." % len(snv_df))
        if len(sv_df) > 0:
            print("Clustering %d SVs..." % len(sv_df))

    Config = ConfigParser.ConfigParser()
    cfg_file = Config.read(cfg)
    dna_gain_class = Config.get('SVclasses', 'dna_gain_class').split(',')
    dna_loss_class = Config.get('SVclasses', 'dna_loss_class').split(',')
    gain_loss_classes = {'dna_gain_class': dna_gain_class, 'dna_loss_class': dna_loss_class}

    cc_sv, cc_snv = None, None
    if len(sv_df) > 0:
        cc_file = '%s/%s_ccube_sv_input.txt' % (out, sample_params['sample'])
        format_svs_for_ccube(sv_df, sample_params, cluster_params, cc_file, gain_loss_classes)

        if os.path.exists(cc_file):
            dirname = os.path.dirname(os.path.abspath(__file__))
            subprocess.call(['Rscript', '%s/cluster_with_ccube.R' % dirname, cc_file,
                             out, sample_params['sample'], '--cores=%d' % threads,
                             '--clusmax=%s' % cluster_params['clus_limit'],
                             '--repeat=%s' % cluster_params['repeat']])

    if len(snv_df) > 0:
        cc_file = '%s/%s_ccube_snv_input.txt' % (out, sample_params['sample'])
        format_snvs_for_ccube(snv_df, sample_params, cluster_params, cc_file)

        if os.path.exists(cc_file):
            dirname = os.path.dirname(os.path.abspath(__file__))
            subprocess.call(['Rscript', '%s/cluster_with_ccube.R' % dirname, cc_file,
                             out, sample_params['sample'], '--cores=%d' % threads,
                             '--clusmax=%s' % cluster_params['clus_limit'],
                             '--repeat=%s' % cluster_params['repeat']])

#    if threads == 1:
#        for run in range(n_runs):
#            snv_run_df = snv_df
#            if subsample > 0 and subsample < len(snv_df):
#                snv_run_df = subsample_snvs(snv_df, subsample, run, ss_seeds, sample, out)
#            cluster_and_process(sv_df,snv_run_df,run,out,sample_params,cluster_params,output_params,seeds)
#        # select the best run based on min BIC
#        if use_map and n_runs > 1:
#            if len(sv_df) > 0:
#                pick_best_run(n_runs, out, sample, ccf_reject, cocluster, fit_metric, cluster_penalty)
#            if len(snv_df) > 0 and not cocluster:
#                pick_best_run(n_runs, out, sample, ccf_reject, cocluster, fit_metric, cluster_penalty, are_snvs=True)
#    else:
#        conc_runs = max(1, n_runs / threads) if n_runs % threads == 0 else (n_runs / threads) + 1
#        for i in range(conc_runs):
#            jobs = []
#            for j in range(threads):
#                run = j + (i * threads)
#                if not run > n_runs-1:
#                    print("Thread %d; cluster run: %d" % (j, run))
#                    snv_run_df = snv_df
#                    if subsample > 0 and subsample < len(snv_df):
#                        snv_run_df = subsample_snvs(snv_df, subsample, run, ss_seeds, sample, out)
#                    process = multiprocessing.Process(target=cluster_and_process,args=(sv_df,snv_run_df,run,out,
#                                                      sample_params,cluster_params,output_params,seeds))
#                    jobs.append(process)
#            for j in jobs:
#                j.start()
#            for j in jobs:
#                j.join()
#
#        while True:
#            if not use_map or n_runs == 1:
#                break
#            if np.all(np.array([not j.is_alive() for j in jobs])):
#                if len(sv_df) > 0:
#                    pick_best_run(n_runs, out, sample, ccf_reject, cocluster, fit_metric, cluster_penalty)
#                if len(snv_df) > 0 and not cocluster:
#                    pick_best_run(n_runs, out, sample, ccf_reject, cocluster, fit_metric, cluster_penalty, are_snvs=True)
#                break
