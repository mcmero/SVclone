from __future__ import print_function

import os
import pandas as pd
import numpy as np
import ConfigParser
import math

from . import cluster
from . import load_data
from . import run_filter as filt
from . import dtypes
from . import write_output
from SVprocess import load_data as svp_load

def get_var_to_assign(var_df, var_filt_df, snvs = False):
    vartype = 'SNVs' if snvs else 'SVs'

    n = int(len(var_df))
    if snvs:
        var_df['support'] = map(float, var_df['var'].values)
    var_df = var_df[var_df.support.values > 0]
    print('Filtered out %d %s with no read support' % (int(n - len(var_df)), vartype))

    if not snvs:
        var_df['gtype1'] = np.array(map(filt.remove_zero_copynumbers, var_df.gtype1.values))
        var_df['gtype2'] = np.array(map(filt.remove_zero_copynumbers, var_df.gtype2.values))

        n = int(len(var_df))
        var_df = var_df[np.logical_or(var_df.gtype1.values!='', var_df.gtype2.values!='')]
        print('Filtered out %d %s with no read support' % (int(n - len(var_df)), vartype))

    var_to_assign = var_df.copy()
    if len(var_filt_df) > 0:
        ids = np.array([])
        if snvs:
            ids = zip(var_df.chrom.values, var_df.pos.map(str).values)
            ids = ['%s:%s' % (x,y) for x,y in ids]
        else:
            ids = zip(var_df.chr1.values, var_df.pos1.values, var_df.dir1,
                      var_df.chr2.values, var_df.pos2.values, var_df.dir2)
            ids = ['%s:%s:%s:%s:%s:%s' % (x,y,z,a,b,c) for x,y,z,a,b,c in ids]

        filt_ids = np.array([])
        if snvs:
            filt_ids = zip(var_filt_df.chrom.values, var_filt_df.pos.values)
            filt_ids = ['%s:%s' % (x,y) for x,y in filt_ids]
        else:
            filt_ids = zip(var_filt_df.chr1.values, var_filt_df.pos1.values, var_filt_df.dir1,
                           var_filt_df.chr2.values, var_filt_df.pos2.values, var_filt_df.dir2)
            filt_ids = ['%s:%s:%s:%s:%s:%s' % (x,y,z,a,b,c) for x,y,z,a,b,c in filt_ids]

        var_to_assign = var_to_assign[[fid not in filt_ids for fid in ids]]

    return var_to_assign

def post_assign_vars(var_df, var_filt_df, rundir, sample, sparams, cparams, snvs = False):
    pa_outdir = '%s_post_assign/' % rundir
    if not os.path.exists(pa_outdir):
        os.makedirs(pa_outdir)

    male     = cparams['male']
    adjusted = cparams['adjusted']
    cnv_pval = cparams['clonal_cnv_pval']
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
    clus_probs = np.array([])

    if len(scs) > 1:
        for i in range(Nvar):
            lls = np.array(np.max(cluster.calc_lik(cn_states[i], sup[i], dep[i], scs.phi.values[0], purity, norm[i])[1]))
            for j in range(1, len(scs)):
                ll = cluster.calc_lik(cn_states[i], sup[i], dep[i], scs.phi.values[j], purity, norm[i])[1]
                lls = np.append(lls, np.max(ll))

            lls_conv = [math.e ** ll for ll in lls]
            ll_probs = [float(lc) / sum(lls_conv) for lc in lls_conv]
            clus_probs = np.concatenate([clus_probs, [ll_probs]], axis=0) if len(clus_probs)>0 else np.array([ll_probs])
            best_clus = np.where(np.max(lls)==lls)[0][0]
            best_clus_list = np.append(best_clus_list, [best_clus])
    else:
        best_clus_list = [0] * Nvar
        clus_probs     = np.array([[1.]] * Nvar)

    best_clus_list = np.array(map(int, best_clus_list))
    phis = scs.phi.values[best_clus_list]

    clusts = scs.clus_id.values
    for clus_idx in best_clus_list:
        clus = clusts[clus_idx]
        idx = scs.index[scs.clus_id.values==clus].values[0]
        scs = scs.set_value(idx, 'size', scs['size'][idx]+1)

    hpd_lo = 'HPD_lo'
    hpd_hi = 'HPD_hi'
    if len(ccert) > 0:
        hpd_lo = ccert.columns.values[-2]
        hpd_hi = ccert.columns.values[-1]
        ccert[hpd_lo] = ccert[hpd_lo]/purity
        ccert[hpd_hi] = ccert[hpd_hi]/purity


    lolim = 0 if snvs else 1
    uplim = 2 if snvs else 7
    ccert_add = var_df[var_df.columns.values[lolim:uplim]]
    ccert_add['most_likely_assignment'] = best_clus_list
    ccert_add['average_ccf'] = np.array(phis)
    ccert_add[hpd_lo] = float('nan')
    ccert_add[hpd_hi] = float('nan')
    ccert = pd.concat([ccert, ccert_add], axis=0)

    # NOTE: clus probabilities are of a different format here
    # CN likelihoods are converted to probabilities,
    # rather than probs being based on MCMC traces
    probs_add = var_df[var_df.columns.values[lolim:uplim]]
    uplim = uplim if snvs else uplim - 1
    prob_cols = probs.columns.values[uplim:]
    for idx,pc in enumerate(prob_cols):
        probs_add[pc] = [round(cp, 2) for cp in clus_probs[:,idx]]

    probs = pd.concat([probs, probs_add], axis=0)
    df = pd.concat([var_filt_df, var_df], axis=0)
    sup = np.append(fsup, sup) if len(var_filt_df) > 0 else sup
    dep = np.append(fdep, dep) if len(var_filt_df) > 0 else dep
    norm = np.append(fnorm, norm) if len(var_filt_df) > 0 else norm
    cn_states = np.append(fcn_states, cn_states) if len(var_filt_df) > 0 else cn_states
    clus_members = np.array([np.where(ccert.most_likely_assignment.values==i)[0] for i in scs.clus_id.values])

    vartype = 'snvs' if snvs else 'svs'
    df.to_csv('%s/../%s_filtered_%s_post_assign.tsv' % (rundir, sample, vartype), sep='\t', index=False, na_rep='')
    print('Writing output to %s' % pa_outdir)
    write_output.write_out_files(df, scs, clus_members, probs, ccert,
                                 pa_outdir, sample, purity, sup, dep,
                                 norm, cn_states, fit, False, cnv_pval, are_snvs = snvs)

def amend_coclus_results(rundir, sample, sparams):
    '''
    correct variant numbers in subclones output file by
    adding post assigned variants SVs and SNVs together
    '''
    purity   = sparams['pi']
    scs, ccert, probs, fit = load_data.get_run_output(rundir, sample, purity)

    pa_outdir = '%s_post_assign/' % rundir
    sv_scs, ccert, probs, fit = load_data.get_run_output(pa_outdir, sample, purity, snvs = False)
    snv_scs, ccert, probs, fit = load_data.get_run_output(pa_outdir, sample, purity, snvs = True)

    correct_sizes = (sv_scs['size'].values - scs['size'].values) + snv_scs['size'].values
    scs['size'] = correct_sizes

    scs['CCF'] = scs.phi.values
    scs['phi'] = scs.phi.values * purity
    scs = scs[['clus_id','size','phi','CCF']]
    rename_cols =  {'clus_id': 'cluster', 'size': 'n_ssms', 'phi': 'proportion'}
    scs = scs.rename(columns = rename_cols)

    scs.to_csv('%s/%s_subclonal_structure.txt' % (pa_outdir, sample), sep='\t', index=False)
    scs.to_csv('%s/snvs/%s_subclonal_structure.txt' % (pa_outdir, sample), sep='\t', index=False)

def run_post_assign(args):

    sample          = args.sample
    cfg             = args.cfg
    out             = sample if args.out == "" else args.out
    snv_file        = args.snv_file
    snv_format      = args.snv_format
    sv_file         = args.sv_file
    cnv_file        = args.cnv_file
    sv_filt_file    = args.sv_filt_file
    snv_filt_file   = args.snv_filt_file
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

    sv_filt_file = '%s/%s_filtered_svs.tsv' % (out, sample) if sv_filt_file == '' else sv_filt_file
    snv_filt_file = '%s/%s_filtered_snvs.tsv' % (out, sample) if snv_filt_file == '' else snv_filt_file

    sv_filt_df  = pd.DataFrame()
    snv_filt_df = pd.DataFrame()

    if os.path.exists(sv_filt_file):
        sv_filt_df = pd.read_csv(sv_filt_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        sv_filt_df = pd.DataFrame(sv_filt_df).fillna('')

    if os.path.exists(snv_filt_file):
        snv_filt_df = pd.read_csv(snv_filt_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        snv_filt_df = pd.DataFrame(snv_filt_df).fillna('')
        snv_filt_df['support'] = map(float, snv_filt_df['var'].values)

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
            snv_df['gtype'] = np.array(map(filt.remove_zero_copynumbers, snv_df.gtype.values))
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
        if len(snv_to_assign) > 0:
            print('Post assigning %d SNVs...' % len(snv_to_assign))
            post_assign_vars(snv_to_assign, snv_filt_df, rundir, sample, sample_params, cluster_params, snvs = True)

    if len(sv_df) > 0 and len(snv_df) > 0:
        amend_coclus_results(rundir, sample, sample_params)

