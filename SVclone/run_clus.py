'''
Facilitates clustering of SVs
'''
from __future__ import print_function

import subprocess
import os
import ConfigParser
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import colorsys
import IPython
import multiprocessing
import time
import shutil

from distutils.dir_util import copy_tree
from collections import OrderedDict
from IPython.core.pylabtools import figsize
from pymc.utils import hpd

from . import cluster
from . import load_data
from . import write_output
from SVprocess import load_data as svp_load

import numpy as np
from numpy import loadtxt
from numpy import zeros

def gen_new_colours(N):
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def plot_clusters(center_trace,clusters,assignments,sup,dep,clus_out_dir):#,alpha_trace):
    fig, axes = plt.subplots(3, 1, sharex=False, sharey=False, figsize=(12.5,8))

    RGB_tuples = gen_new_colours(len(clusters))

    axes[0].set_ylim([0,2])
    axes[0].set_title("Trace of $\phi_k$")

    for idx,clus in enumerate(clusters):
        axes[0].plot(center_trace[:, clus], label="trace of center %d" % clus, c=RGB_tuples[idx], lw=1)
        # print('phi_%d ~ %f' % (idx,np.mean(center_trace[2000:,clus])))

    leg = axes[0].legend(loc="upper right")
    leg.get_frame().set_alpha(0.7)
    axes[1].set_title("Raw VAFs")
#    axes[2].set_title("Alpha trace")
#    axes[2].plot(alpha_trace, label="alpha trace", lw=1)
    
    for idx,clus in enumerate(clusters):
        clus_idx = np.array(assignments)==clus
        sup_clus = sup[clus_idx]
        dep_clus = dep[clus_idx]
        axes[1].hist(sup_clus/dep_clus,bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])
        #axes[2].hist(((s1/(n1+s1)*pl)/pi),bins=np.array(range(0,100,2))/100.,alpha=0.75,color=RGB_tuples[idx])

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

def merge_clusters(clus_out_dir,clus_info,clus_merged,clus_members,merged_ids,sup,dep,cn_states,sparams,cparams,subclone_diff,clus_limit,phi_limit):
    if len(clus_info)==1:
        return (clus_info,clus_members)

    clus_info = clus_info.sort('phi',ascending=False).copy()
    clus_info.index = range(0,len(clus_info))
    to_del = []
    for idx,ci in clus_info.iterrows():
        if idx+1 >= len(clus_info):
            clus_merged.loc[idx] = clus_info.loc[idx]
            break
        cn = clus_info.loc[idx+1]
        if abs(ci.phi - float(cn.phi)) < subclone_diff:
            print("\nReclustering similar clusters...")
            new_size = ci['size'] + cn['size']

            new_members = np.concatenate([clus_members[idx],clus_members[idx+1]])
            mcmc, map_ = cluster.cluster(sup[new_members],dep[new_members],cn_states[new_members],len(new_members),sparams,cparams,clus_limit,phi_limit)
            trace = mcmc.trace("phi_k")[:]

            phis = mean_confidence_interval(trace,cparams['hpd_alpha'])
            clus_merged.loc[idx] = np.array([ci.clus_id,new_size,phis[0],phis[1],phis[2]])            
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
        clus_merged = clus_merged.sort('phi',ascending=False)
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
        return merge_clusters(clus_out_dir,clus_merged,new_df,clus_members,merged_ids,sup,dep,cn_states,sparams,cparams,subclone_diff,clus_limit,phi_limit)

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

def post_process_clusters(mcmc,sv_df,snv_df,clus_out_dir,sup,dep,cn_states,sparams,cparams,output_params,map_):

    merge_clusts  = cparams['merge_clusts']
    subclone_diff = cparams['subclone_diff']
    phi_limit     = cparams['phi_limit']
    merge_clusts  = cparams['merge_clusts']
    smc_het       = output_params['smc_het']
    write_matrix  = output_params['write_matrix']
    plot          = output_params['plot']

    # fit to data frame
    run_fit = pd.DataFrame()
    if map_ is not None:
        run_fit = pd.DataFrame([['BIC', map_.BIC], ['AIC', map_.AIC]])

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
    
    if len(clus_info) < 1:
        print("Warning! Could not converge on any major SV clusters. Skipping.\n")
        return None

    # get cluster means
    center_trace = mcmc.trace("phi_k")[:]
    
    phis = np.array([mean_confidence_interval(center_trace[:,cid],cparams['hpd_alpha']) for cid in clus_info.clus_id.values])
    clus_info['phi'] = phis[:,0]
    clus_info['95p_HPD_lo'] = phis[:,1]
    clus_info['95p_HPD_hi'] = phis[:,2]
    
    clus_ids = clus_info.clus_id.values
    clus_members = np.array([np.where(np.array(clus_max_prob)==i)[0] for i in clus_ids])

    col_names = map(lambda x: 'cluster'+str(x),clus_ids)
    df_probs = pd.DataFrame(clus_counts,dtype=float)[clus_ids].fillna(0)
    df_probs = df_probs.apply(lambda x: x/sum(x),axis=1)
    df_probs.columns = col_names

    # cluster certainty
    clus_max_df = pd.DataFrame(clus_max_prob,columns=['most_likely_assignment'])
    phi_cols = ["average_ccf","95p_HPD_lo","95p_HPD_hi"]
    phi_matrix = pd.DataFrame(phis[:],index=clus_ids,columns=phi_cols).loc[clus_max_prob]
    phi_matrix.index = range(len(phi_matrix))
    ccert = clus_max_df.join(phi_matrix)
    clus_info.index = range(len(clus_info))

    print('\n\n')
    print(clus_info)
    print('Compiling and writing output...')

    dump_out_dir = clus_out_dir
    if len(snv_df)>0 and len(sv_df)==0:
        # snvs only trace output
        dump_out_dir = '%s/snvs'%clus_out_dir        
    trace_out = 'premerge_' if merge_clusts else ''
    trace_out = '%s/%s'%(dump_out_dir,trace_out)
    write_output.dump_trace(clus_info,center_trace,trace_out+'phi_trace.txt')
    write_output.dump_trace(clus_info,z_trace,trace_out+'z_trace.txt')
    
    # cluster plotting
    if plot:
        alpha_trace = mcmc.trace("alpha")[:]
        plot_clusters(center_trace,clus_idx,clus_max_prob,sup,dep,clus_out_dir)
    
    # merge clusters
    if len(clus_info)>1 and merge_clusts:        
        clus_merged = pd.DataFrame(columns=clus_info.columns,index=clus_info.index)
        clus_merged, clus_members, merged_ids  = merge_clusters(clus_out_dir,clus_info,clus_merged,\
                clus_members,[],sup,dep,cn_states,sparams,cparams,subclone_diff,clus_limit,phi_limit)
        
        if len(clus_merged)!=len(clus_info):
            clus_info = clus_merged
            df_probs, ccert = merge_results(clus_merged, merged_ids, df_probs, ccert)    

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

        snv_sup  = sup[:len(snv_df)]
        snv_dep  = dep[:len(snv_df)]
        snv_cn_states = cn_states[:len(snv_df)]
        write_output.write_out_files(snv_df,clus_info,snv_members,
                snv_probs,snv_ccert,clus_out_dir,sparams['sample'],sparams['pi'],snv_sup,
                snv_dep,snv_cn_states,run_fit,z_trace,smc_het,write_matrix,are_snvs=True)
    
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

        sv_sup  = sup[lb:lb+len(sv_df)]
        sv_dep  = dep[lb:lb+len(sv_df)]
        sv_cn_states = cn_states[lb:lb+len(sv_df)]
        write_output.write_out_files(sv_df,clus_info,sv_members,
                    sv_probs,sv_ccert,clus_out_dir,sparams['sample'],sparams['pi'],sv_sup,
                    sv_dep,sv_cn_states,run_fit,z_trace,smc_het,write_matrix)

def cluster_and_process(sv_df, snv_df, run, out_dir, sample_params, cluster_params, output_params):
    # set random seed
    seed = int(round((time.time() - int(time.time())) * 10000 * (run+1)))
    np.random.seed(seed)

    clus_out_dir = '%s/run%d'%(out_dir, run)
    if not os.path.exists(clus_out_dir):
        os.makedirs(clus_out_dir)

    Nvar, sup, dep, cn_states = None, None, None, None
    if cluster_params['cocluster'] and len(sv_df)>0 and len(snv_df)>0:
        # coclustering
        sup, dep, cn_states, Nvar = load_data.get_snv_vals(snv_df)
        sv_sup, sv_dep, sv_cn_states, sv_Nvar = load_data.get_sv_vals(sv_df, cluster_params['adjusted'])
        sup = np.append(sup, sv_sup)
        dep = np.append(dep, sv_dep)
        cn_states = pd.concat([pd.DataFrame(cn_states),pd.DataFrame(sv_cn_states)])[0].values
        Nvar = Nvar + sv_Nvar
        print("Coclustering %d SVs & %d SNVs..." % (len(sv_df), len(snv_df)))
        mcmc, map_ = cluster.cluster(sup, dep, cn_states, Nvar, sample_params,
                                     cluster_params, cluster_params['phi_limit'])
        post_process_clusters(mcmc, sv_df, snv_df, clus_out_dir, sup, dep, cn_states,
                          sample_params, cluster_params, output_params, map_)

    elif len(sv_df) > 0 or len(snv_df) > 0:
        # no coclustering
        if len(snv_df) > 0:
            if not os.path.exists('%s/snvs'%clus_out_dir):
                os.makedirs('%s/snvs'%clus_out_dir)
            sup,dep,cn_states,Nvar = load_data.get_snv_vals(snv_df)
            print('Clustering %d SNVs...' % len(snv_df))
            mcmc, map_ = cluster.cluster(sup, dep, cn_states, Nvar, sample_params,
                                         cluster_params, cluster_params['phi_limit'])
            post_process_clusters(mcmc, pd.DataFrame(), snv_df, clus_out_dir, sup, dep, cn_states,
                              sample_params, cluster_params, output_params, map_)
        if len(sv_df) > 0:
            sup, dep, cn_states, Nvar = load_data.get_sv_vals(sv_df,cluster_params['adjusted'])
            print('Clustering %d SVs...' % len(sv_df))
            mcmc, map_ = cluster.cluster(sup, dep, cn_states, Nvar, sample_params,
                                         cluster_params, cluster_params['phi_limit'])
            post_process_clusters(mcmc, sv_df, pd.DataFrame(), clus_out_dir, sup, dep, cn_states,
                              sample_params, cluster_params, output_params, map_)

    else:
        raise ValueError("No valid variants to cluster!")

def string_to_bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def run_clustering(args):
    
    snv_file        = args.snv_file
    sv_file         = args.sv_file
    sample          = args.sample
    out             = args.outdir
    param_file      = args.param_file
    pp_file         = args.pp_file
    cfg             = args.cfg

    out = sample if out == "" else out
    Config = ConfigParser.ConfigParser()
    cfg_file = Config.read(cfg)

    if len(cfg_file)==0:
        raise ValueError('No configuration file found')

    shape  = Config.get('BetaParameters', 'alpha')
    scale  = Config.get('BetaParameters', 'beta')
    beta   = ','.join([str(shape),str(scale)])

    phi_limit       = float(Config.get('ClusterParameters', 'phi_limit'))
    clus_limit      = int(Config.get('ClusterParameters', 'clus_limit'))
    subclone_diff   = float(Config.get('ClusterParameters', 'subclone_diff'))
    hpd_alpha       = float(Config.get('ClusterParameters', 'hpd_alpha'))

    n_runs          = int(Config.get('ClusterParameters', 'n_runs'))
    n_iter          = int(Config.get('ClusterParameters', 'n_iter'))
    burn            = int(Config.get('ClusterParameters', 'burn'))
    thin            = int(Config.get('ClusterParameters', 'thin'))
    threads         = int(Config.get('ClusterParameters', 'threads'))

    use_map         = string_to_bool(Config.get('ClusterParameters', 'map'))
    merge_clusts    = string_to_bool(Config.get('ClusterParameters', 'merge'))
    cocluster       = string_to_bool(Config.get('ClusterParameters', 'cocluster'))
    adjusted       = string_to_bool(Config.get('ClusterParameters', 'adjusted'))

    plot            = string_to_bool(Config.get('OutputParameters', 'plot'))
    smc_het         = string_to_bool(Config.get('OutputParameters', 'smc_het'))
    write_matrix    = string_to_bool(Config.get('OutputParameters', 'coclus_matrix'))

    if out!='' and not os.path.exists(out):
        os.makedirs(out)

    pi, pl = svp_load.get_purity_ploidy(pp_file, sample, out)
    rlen, insert, std = svp_load.get_read_params(param_file, sample, out)

    sample_params  = { 'sample': sample, 'ploidy': pl, 'pi': pi, 'rlen': rlen, 'insert': insert }
    cluster_params = { 'n_iter': n_iter, 'burn': burn, 'thin': thin, 'beta': beta, 'use_map': use_map, 'hpd_alpha': hpd_alpha,
                       'merge_clusts': merge_clusts, 'adjusted': adjusted, 'phi_limit': phi_limit, 'clus_limit': clus_limit,
                       'subclone_diff': subclone_diff, 'cocluster': cocluster }
    output_params  = { 'plot': plot, 'write_matrix': write_matrix, 'smc_het': smc_het }
    
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

    clus_info,center_trace,ztrace,clus_members = None,None,None,None

    if threads == 1:
        for run in range(n_runs):
            cluster_and_process(sv_df,snv_df,run,out,sample_params,cluster_params,output_params)
    else:
        conc_runs = max(1, n_runs / threads) if n_runs % threads == 0 else (n_runs / threads) + 1
        for i in range(conc_runs):
            jobs = []

            for j in range(threads):

                run = j + (i * threads)
                if not run > n_runs-1:
                    print("Thread %d; cluster run: %d" % (j, run))
                    process = multiprocessing.Process(target=cluster_and_process,args=(sv_df,snv_df,run,out,
                                                      sample_params,cluster_params,output_params))
                    jobs.append(process)

            for j in jobs:
                j.start()
            for j in jobs:
                j.join()

    # select the best run based on min BIC
    if use_map:
        bics = []
        if len(sv_df) > 0:
            for run in range(n_runs):
                fit_file = '%s/run%d/%s_fit.txt' % (out, run, sample)
                fit = pd.read_csv(fit_file, delimiter='\t', dtype=None, header=None)
                bics.append(fit.loc[0][1])
            bics = np.array(bics)
            min_bic = np.where(min(bics) == bics)[0][0]

            print('Selecting run %d as best run for SVs' % min_bic)
            best_run = '%s/run%d' % (out, min_bic)
            copy_tree(best_run, '%s/best_run_svs' % out)
            snv_folder = '%s/best_run_svs/snvs' % out
            if os.path.exists(snv_folder):
                shutil.rmtree(snv_folder)

        bics = []
        if len(snv_df) > 0:
            for run in range(n_runs):
                fit_file = '%s/run%d/snvs/%s_fit.txt' % (out, run, sample)
                fit = pd.read_csv(fit_file, delimiter='\t', dtype=None, header=None)
                bics.append(fit.loc[0][1])
            bics = np.array(bics)
            min_bic = np.where(min(bics) == bics)[0][0]

            print('Selecting run %d as best run for SNVs' % min_bic)
            best_run = '%s/run%d/snvs' % (out, min_bic)
            best_run_dest = '%s/best_run_snvs' % out
            copy_tree(best_run, best_run_dest)

            if cocluster:
                shutil.copyfile('%s/run%d/phi_trace.txt' % (out, min_bic), '%s/phi_trace.txt' % best_run_dest)
                shutil.copyfile('%s/run%d/z_trace.txt' % (out, min_bic), '%s/z_trace.txt' % best_run_dest)
            else:
                shutil.copyfile('%s/run%d/snvs/phi_trace.txt' % (out, min_bic), '%s/phi_trace.txt' % best_run_dest)
                shutil.copyfile('%s/run%d/snvs/z_trace.txt' % (out, min_bic), '%s/z_trace.txt' % best_run_dest)

            shutil.copyfile('%s/run%d/cluster_trace_hist.png' % (out, min_bic), 
                            '%s/cluster_trace_hist.png' % best_run_dest)
