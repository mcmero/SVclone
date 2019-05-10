'''
Facilitates clustering of SVs
'''
from __future__ import print_function

import subprocess
import os
import pandas as pd
import numpy as np
import random
import ConfigParser

from . import load_data
from SVprocess import svp_load_data as svp_load


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
                         'major_cn_sub1': cn_states[0],
                         'minor_cn_sub1': cn_states[1],
                         'total_cn_sub1': cn_states[0] + cn_states[1],
                         'frac_cn_sub1': cn_states[2],
                         'major_cn_sub2': cn_states[3],
                         'minor_cn_sub2': cn_states[4],
                         'total_cn_sub2': cn_states[3] + cn_states[4],
                         'frac_cn_sub2': cn_states[5]})

    cc_df.to_csv(cc_file, sep='\t', index=False)

def format_svs_for_ccube(df, sparams, cparams, cc_file):
    '''
    prepare ccube input for SVs
    '''
    dl = len(df)
    df = df[np.logical_and(df.gtype1!='', df.gtype2!='')]
    if len(df) != dl:
        print('Filtered out %d SVs with missing copy-number states' % (dl - len(df)))

    sup, norm1, norm2, Nvar, norm_cn = load_data.get_sv_vals(df, cparams)
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
                         'ref_counts1': norm1,
                         'total_counts1': sup + norm1,
                         'vaf1': sup / (sup + norm1),
                         'major_cn1_sub1': cn_states1[0],
                         'minor_cn1_sub1': cn_states1[1],
                         'total_cn1_sub1': cn_states1[0] + cn_states1[1],
                         'frac_cn1_sub1': cn_states1[2],
                         'major_cn1_sub2': cn_states1[3],
                         'minor_cn1_sub2': cn_states1[4],
                         'total_cn1_sub2': cn_states1[3] + cn_states1[4],
                         'frac_cn1_sub2': cn_states1[5],
                         'var_counts2': sup,
                         'ref_counts2': norm2,
                         'total_counts2': sup + norm2,
                         'vaf2': sup / (sup + norm2),
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
    XX              = args.XX
    XY              = args.XY
    subsample       = args.subsample
    ss_seeds        = args.ss_seeds

    if out!='' and not os.path.exists(out):
        os.makedirs(out)

    sample_params, cluster_params = \
                   load_data.get_params_cluster_step(sample, cfg, out, pp_file, param_file, XX, XY)
    ss_seeds = get_seeds(ss_seeds, 1)

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

    sv_to_sim = cluster_params['sv_to_sim']
    if sv_to_sim > 0:
        print('Simulating %d SV...' % (len(sv_df) * sv_to_sim))
        nvar = len(sv_df)
        for i in range(nvar):
            for j in range(sv_to_sim):
                sv_df.loc[len(sv_df)] = simu_sv(sv_df.loc[i])

    threads = cluster_params['threads']
    cc_sv, cc_snv = None, None
    if len(sv_df) > 0:
        cc_file = '%s/%s_ccube_sv_input.txt' % (out, sample_params['sample'])
        format_svs_for_ccube(sv_df, sample_params, cluster_params, cc_file)

        if os.path.exists(cc_file):
            dirname = os.path.dirname(os.path.abspath(__file__))
            subprocess.call(['Rscript', '%s/cluster_with_ccube.R' % dirname, cc_file,
                             out, sample_params['sample'], '--cores=%d' % threads,
                             '--clusmax=%s' % cluster_params['clus_limit'],
                             '--repeat=%s' % cluster_params['repeat'],
                             '--maxiter=%s' % cluster_params['n_iter']])

    if len(snv_df) > 0:
        cc_file = '%s/%s_ccube_snv_input.txt' % (out, sample_params['sample'])
        format_snvs_for_ccube(snv_df, sample_params, cluster_params, cc_file)

        if os.path.exists(cc_file):
            dirname = os.path.dirname(os.path.abspath(__file__))
            subprocess.call(['Rscript', '%s/cluster_with_ccube.R' % dirname, cc_file,
                             out, sample_params['sample'], '--cores=%d' % threads,
                             '--clusmax=%s' % cluster_params['clus_limit'],
                             '--repeat=%s' % cluster_params['repeat'],
                             '--maxiter=%s' % cluster_paramx['n_iter']])
