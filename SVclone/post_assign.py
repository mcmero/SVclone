from __future__ import print_function

import os
import pandas as pd
import numpy as np
import ConfigParser
import ipdb

from . import cluster
from . import load_data
from . import run_filter as filt
from . import dtypes
from . import write_output
from SVprocess import load_data as svp_load

def get_var_to_assign(var_df, var_filt_df, snvs = False):
    vartype = 'SNVs' if snvs else 'SVs'
    var_to_assign = var_df

    if len(var_filt_df) > 0:
        ids = np.array([])
        if snvs:
            ids = zip(snv_df.chrom.values, snv_df.pos.values)
        else:
            ids = zip(var_df.chr1.values, var_df.pos1.values, var_df.dir1,
                      var_df.chr2.values, var_df.pos2.values, var_df.dir2)
        ids = [':'.join([str(y) for y in x]) for x in ids]

        filt_ids = np.array([])
        if snvs:
            filt_ids = zip(snv_filt_df.chrom.values, snv_filt_df.pos.values),
        else:
            filt_ids = zip(var_filt_df.chr1.values, var_filt_df.pos1.values, var_filt_df.dir1,
                           var_filt_df.chr2.values, var_filt_df.pos2.values, var_filt_df.dir2)
        filt_ids = [':'.join([str(y) for y in x]) for x in filt_ids]

        var_to_assign = var_df[[id not in filt_ids for id in ids]]

    n = int(len(var_to_assign))
    if snvs:
        var_to_assign['support'] = map(int, var_to_assign['var'].values)
    var_to_assign = var_to_assign[var_to_assign.support.values > 0]
    print('Filtered out %d %s with no read support' % (int(n - len(var_to_assign)), vartype))

    return var_to_assign

def post_assign_vars(var_df, var_filt_df, rundir, sample, sparams, cparams, snvs = False):
    pa_dir = '%s_post_assign' % rundir if not snvs else '%s_post_assign/snvs' % rundir
    if not os.path.exists(pa_dir):
        os.makedirs(pa_dir)
    print('Writing output to %s' % pa_dir)

    male     = cparams['male']
    adjusted = cparams['adjusted']
    purity   = sparams['pi']
    ploidy   = sparams['ploidy']
    mean_cov = sparams['mean_cov']

    scs, ccert, probs, fit = load_data.get_run_output(rundir, sample, purity, snvs)
    sup,   dep,  cn_states,  Nvar,  norm = None, None, None, None, None
    fsup, fdep, fcn_states, fNvar, fnorm = None, None, None, None, None

    if snvs:
        sup, dep, cn_states, Nvar, norm = load_data.get_snv_vals(var_df, male)
        if len(var_filt_df) > 0:
            fsup, fdep, fcn_states, fNvar, fnorm = load_data.get_snv_vals(var_filt_df, male)
    else:
        sup, dep, cn_states, Nvar, norm = load_data.get_sv_vals(var_df, adjusted, male)
        if len(var_filt_df) > 0:
            fsup, fdep, fcn_states, fNvar, fnorm = load_data.get_sv_vals(var_filt_df, adjusted, male)

    best_clus_list = np.array([])
    if len(scs) > 0:
        for i in range(Nvar):
            lls = np.array(np.max(cluster.calc_lik(cn_states[i], sup[i], dep[i], scs.phi.values[0], purity, norm[i])[1]))
            for j in range(1, len(scs)):
                ll = cluster.calc_lik(cn_states[i], sup[i], dep[i], scs.phi.values[j], purity, norm[i])[1]
                lls = np.append(lls, np.max(ll))

            best_clus = np.where(np.max(lls)==lls)[0][0]
            best_clus_list = np.append(best_clus_list, [best_clus])
    else:
        best_clust_list = [0] * Nvar

    best_clus_list = np.array(map(int, best_clus_list))
    phis = scs.phi.values[best_clus_list]

    for clus in best_clus_list:
        idx = scs.index[scs.clus_id.values==clus].values[0]
        scs = scs.set_value(idx, 'size', scs['size'][idx]+1)

    lim = 3 if snvs else 6
    ccert_add = var_df[var_df.columns.values[:lim]]
    ccert_add['most_likely_assignment'] = best_clus_list
    ccert_add['average_ccf'] = np.array(phis)
    ccert_add['95_HPD_lo'] = float('nan')
    ccert_add['95_HPD_hi'] = float('nan')
    ccert = pd.concat([ccert, ccert_add])

    probs_add = var_df[var_df.columns.values[:lim]]
    prob_cols = probs.columns.values[lim:]
    for i in prob_cols:
        probs_add[i] = float('nan')

    probs = pd.concat([probs, probs_add])
    df = pd.concat([var_filt_df, var_df])
    sup = np.append(fsup, sup) if len(var_filt_df) > 0 else sup
    dep = np.append(fdep, dep) if len(var_filt_df) > 0 else dep
    norm = np.append(fnorm, norm) if len(var_filt_df) > 0 else norm
    cn_states = np.append(fcn_states, cn_states) if len(var_filt_df) > 0 else cn_states
    clus_members = ccert.most_likely_assignment.values

    # CNV pval cutoff set to 0 (which means always select based on phi, not clonal)
    write_output.write_out_files(df, scs, clus_members, probs, ccert,
                                 pa_dir, sample, purity, sup, dep,
                                 norm, cn_states, fit, False, 0, are_snvs = snvs)

def run_post_assign(args):

    sample          = args.sample
    cfg             = args.cfg
    out             = sample if args.out == "" else args.out
    snv_file        = args.snv_file
    snv_format      = args.snv_format
    sv_file         = args.sv_file
    cnv_file        = args.cnv_file
    run             = args.run
    gml             = args.germline
    XX              = args.XX
    XY              = args.XY

    if out == '':
        out = sample

    Config = ConfigParser.ConfigParser()
    cfg_file = Config.read(cfg)
    if len(cfg_file)==0:
        raise ValueError('No configuration file found')

    sv_offset = int(Config.get('FilterParameters', 'sv_offset'))
    gl_th     = int(Config.get('FilterParameters', 'germline_threshold'))

    if out != '' and not os.path.exists(out):
        raise ValueError('Specified output directory does not exist!')

    if run == 'best':
        run_dirs = next(os.walk(out))[1]
        run_dirs = [run for run in run_dirs if run.startswith('best_run') and not run.endswith('post_assign')]
        if len(run_dirs) > 1:
            raise ValueError('More than 1 best run directory exists! Please specify which run to use.')
        elif len(run_dirs) < 1:
            raise ValueError('No best run directories exist! Please specify which run to use.')
        run = run_dirs[0]

    rundir = '%s/%s' % (out, run)
    if not os.path.exists(rundir):
        raise OSError('Specified run directory does not exist!')

    if snv_file == '' and sv_file == '':
        raise ValueError('No variants specified!')

    if snv_file != '' and not os.path.exists(snv_file):
        raise OSError('Specified SNV file does not exist!')

    if sv_file != '' and not os.path.exists(sv_file):
        raise OSError('Specified SV file does not exist!')

    sv_filt_file = '%s/%s_filtered_svs.tsv' % (out, sample)
    snv_filt_file = '%s/%s_filtered_snvs.tsv' % (out, sample)

    sv_filt_df  = pd.DataFrame()
    snv_filt_df = pd.DataFrame()

    if os.path.exists(sv_filt_file):
        sv_filt_df = pd.read_csv(sv_filt_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        sv_filt_df = pd.DataFrame(sv_filt_df).fillna('')

    if os.path.exists(snv_filt_file):
        snv_filt_df = pd.read_csnv(snv_filt_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        snv_filt_df = pd.DataFrame(snv_filt_df).fillna('')

    if len(sv_filt_df) == 0 and len(snv_filt_df) == 0:
        raise ValueError('Output directory filtered variant files do not exist or are empty!')

    param_file = '%s/read_params.txt' % out
    pp_file = '%s/purity_ploidy.txt' % out

    sample_params, cluster_params, output_params = \
                   load_data.get_params_cluster_step(sample, cfg, out, pp_file, param_file, XX, XY)
    rlen, insert, insert_std = svp_load.get_read_params(param_file, sample, out)
    purity = sample_params['pi']
    ploidy = sample_params['ploidy']

    sv_df  = pd.DataFrame()
    snv_df = pd.DataFrame()

    if snv_file != '':
        if snv_format == 'sanger':
            snv_df = load_data.load_snvs_sanger(snv_file)
        elif snv_format == 'mutect':
            snv_df = load_data.load_snvs_mutect(snv_file, sample)
        elif snv_format == 'mutect_callstats':
            snv_df = load_data.load_snvs_mutect_callstats(snv_file)
        elif snv_format == 'consensus':
            snv_df = load_data.load_snvs_consensus(snv_file)

    if sv_file != "":
        sv_df = load_data.load_svs(sv_file)
        if gml!="":
            sv_df = filt.filter_germline(gml, sv_df, rlen, insert, gl_th)

    if cnv_file != "":
        cnv_df = load_data.load_cnvs(cnv_file)
        strict_cnv_filt = False

        if len(sv_df)>0:
            print('Matching copy-numbers for SVs...')
            sv_df = filt.match_copy_numbers(sv_df,cnv_df,strict_cnv_filt,sv_offset)
            sv_df = filt.match_copy_numbers(sv_df,cnv_df,strict_cnv_filt,sv_offset,\
                    ['chr2','pos2','dir2','classification','pos1'],'gtype2')

        if len(snv_df)>0:
            print('Matching copy-numbers for SNVs...')
            snv_df = filt.match_snv_copy_numbers(snv_df,cnv_df,strict_cnv_filt)
            n = len(snv_df)
            snv_df = snv_df[snv_df.gtype != '']
            print('Filtered out %d SNVs with no copy-numbers' % (n-len(snv_df)))
    else:
        print('No CNV input defined, assuming all loci major/minor allele copy-numbers are ploidy/2')

        maj_allele = round(ploidy/2) if round(ploidy) > 1 else 1
        min_allele = round(ploidy/2) if round(ploidy) > 1 else 0
        default_gtype = '%d,%d,1.0' % (maj_allele, min_allele)

        if len(sv_df)>0:
            sv_df['gtype1'] = default_gtype
            sv_df['gtype2'] = default_gtype

        if len(snv_df)>0:
            snv_df['gtype'] = default_gtype

    if len(sv_df) > 0:
        sv_to_assign = get_var_to_assign(sv_df, sv_filt_df)
        if len(sv_to_assign) > 0:
            sv_to_assign = filt.adjust_sv_read_counts(sv_to_assign, purity, ploidy, 0, rlen, Config)
            print('Post assigning %d SVs...' % len(sv_to_assign))
            post_assign_vars(sv_to_assign, sv_filt_df, rundir, sample, sample_params, cluster_params)

    if len(snv_df) > 0:
        snv_to_assign = get_var_to_assign(snv_df, snv_filt_df, snvs = True)
        print('Post assigning %d SNVs...' % len(snv_to_assign))
        post_assign_vars(snv_to_assign, snv_filt_df, rundir, sample, sample_params, cluster_params, snvs = True)

