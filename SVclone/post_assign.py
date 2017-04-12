from __future__ import print_function

import os
import pandas as pd
import numpy as np
import ConfigParser
import math

from shutil import copytree, ignore_patterns, rmtree

from . import cluster
from . import run_clus
from . import load_data
from . import run_filter as filt
from . import dtypes
from . import write_output
from SVprocess import load_data as svp_load

def get_var_ids(var_df, snvs):
    var_ids = np.array([])

    if snvs:
        var_ids = zip(var_df.chrom.values, var_df.pos.values)
        var_ids = ['%s:%s' % (x,y) for x,y in var_ids]
    else:
        var_ids = zip(var_df.chr1.values, var_df.pos1.values, var_df.dir1,
                  var_df.chr2.values, var_df.pos2.values, var_df.dir2)
        var_ids = ['%s:%s:%s:%s:%s:%s' % (x,y,z,a,b,c) for x,y,z,a,b,c in var_ids]

    return(var_ids)

def get_var_to_assign(var_df, var_filt_df, snvs = False):
    vartype = 'SNVs' if snvs else 'SVs'

    n = int(len(var_df))
    if snvs:
        var_df['support'] = map(float, var_df['var'].values)
    var_df = var_df[var_df.support.values > 0]
    n_filt_out = int(n - len(var_df))
    if n_filt_out > 0:
        print('Filtered out %d %s with no read support' % (n_filt_out, vartype))

    if not snvs:
        var_df['gtype1'] = np.array(map(filt.remove_zero_copynumbers, var_df.gtype1.values))
        var_df['gtype2'] = np.array(map(filt.remove_zero_copynumbers, var_df.gtype2.values))

        n = int(len(var_df))
        var_df = var_df[np.logical_or(var_df.gtype1.values!='', var_df.gtype2.values!='')]
        n_filt_out = int(n - len(var_df))
        if n_filt_out > 0:
            print('Filtered out %d %s with missing copy-number states' % (int(n - len(var_df)), vartype))

    var_to_assign = var_df.copy()
    if len(var_filt_df) > 0:
        ids = get_var_ids(var_df, snvs)
        filt_ids = get_var_ids(var_filt_df, snvs)

        var_to_assign = var_to_assign[[fid not in filt_ids for fid in ids]]

    return var_to_assign

def get_ll_probs(s, d, n, cn, purity, scs):
    lls = np.array([np.max(cluster.calc_lik(cn, s, d, scs.phi.values[0], purity, n)[1])])
    for j in range(1, len(scs)):
        ll = cluster.calc_lik(cn, s, d, scs.phi.values[j], purity, n)[1]
        lls = np.append(lls, np.max(ll))
    if len(lls) > 1:
        lls_conv = [math.e ** ll for ll in lls]
        ll_probs = [float(lc) / sum(lls_conv) for lc in lls_conv]
    else:
        ll_probs = [1.]

    return lls, ll_probs

def fix_subsampling_discrepancy(var_df, var_filt_df, filt_ids, ccert_ids):
    n_to_assign = len(var_df)
    var_df = pd.concat([var_df, var_filt_df])
    var_in_ccert = np.array([var_id in ccert_ids for var_id in filt_ids])

    var_filt_df = var_filt_df[var_in_ccert]
    to_assign = np.concatenate([ np.array([True] * n_to_assign, dtype=bool),
                                 np.invert(var_in_ccert) ])

    var_df.index = range(len(var_df))
    var_df = var_df[to_assign]
    return(var_df, var_filt_df)

def filter_clusters(scs, clus_th, Nvar, sup, dep, norm, cn_states, purity):
    '''
    filter out clusters of low proportion
    '''

    best_clus_list = np.array([])
    clus_probs = np.array([])

    percent_filt = (scs['size']/sum(scs['size'])).values > clus_th['percent']
    abs_filt     = scs['size'].values > clus_th['absolute']
    orig_scs     = scs.copy()
    scs          = scs[np.logical_and(percent_filt, abs_filt)]

    if len(scs) == 0:
        scs = orig_scs
        print("WARNING: all clusters fail minimum size requirements! No cluster filtering has been performed.")

    if len(scs) > 1 and Nvar > 0:
        for i in range(Nvar):
            lls, ll_probs = get_ll_probs(sup[i], dep[i], norm[i], cn_states[i], purity, scs)
            clus_probs = np.concatenate([clus_probs, [ll_probs]], axis=0) if len(clus_probs)>0 else np.array([ll_probs])
            best_clus = np.where(np.max(lls)==lls)[0][0]
            best_clus_list = np.append(best_clus_list, [best_clus])
    elif Nvar > 0:
        best_clus_list = [0] * Nvar
        clus_probs     = np.array([[1.]] * Nvar)

    best_clus_list = np.array(map(int, best_clus_list))
    phis = scs.phi.values[best_clus_list] if Nvar > 0 else None

    return scs, orig_scs, best_clus_list, clus_probs, phis

def reclassify_filtered_out_vars(scs, orig_scs, probs, ccert, fsup, fdep, fnorm, fcn_states, purity):
    '''
    reclassify minor cluster variants if filtered out
    '''
    hpd_lo = ccert.columns.values[-2]
    hpd_hi = ccert.columns.values[-1]

    clus_fo = np.array([clus for clus in orig_scs.clus_id.values if clus not in scs.clus_id.values])
    # drops cols from probabilities dataframe
    for i in clus_fo:
        probs = probs.drop('cluster%d' % i, 1)

    mla = ccert.most_likely_assignment.values
    which_reclass = np.where(np.array([cv in clus_fo for cv in mla]))[0]
    if len(which_reclass) > 0:
        not_reclass = np.array([ccidx not in which_reclass for ccidx in ccert.index.values])
        print('Reassigning %d variants from filtered out clusters' % len(which_reclass))
        for i in which_reclass:
            # find most likely cluster for variant
            lls, ll_probs = get_ll_probs(fsup[i], fdep[i], fnorm[i], fcn_states[i], purity, scs)

            # pos refers to index position, while array_idx is the array's index (this may differ)
            best_clus_pos = np.where(np.max(lls)==lls)[0][0]
            best_clus_array_idx = scs.index[scs.clus_id.values==scs.clus_id.values[best_clus_pos]].values[0]

            # change variant's most likely assignment to new membership
            new_clus = scs.clus_id.values[best_clus_pos]
            ccert = ccert.set_value(i, 'most_likely_assignment', new_clus)

            # update cluster proportion
            ccert = ccert.set_value(i, 'average_ccf', scs.phi.values[best_clus_pos])
            tmp = ccert.loc[not_reclass]
            lo_val = tmp[tmp.most_likely_assignment==new_clus][hpd_lo].values[0]
            hi_val = tmp[tmp.most_likely_assignment==new_clus][hpd_hi].values[0]
            ccert = ccert.set_value(i, hpd_lo, lo_val)
            ccert = ccert.set_value(i, hpd_hi, hi_val)

            # update cluster assignment probability
            pa_clus_cols = ['cluster%d' % cid for cid in scs.clus_id.values]
            if len(pa_clus_cols) > 1:
                for idx,pac in enumerate(pa_clus_cols):
                    probs = probs.set_value(i, pac, round(ll_probs[idx], 2))
            else:
                probs.set_value(i, pa_clus_cols[0], 1.)

            # increment subclonal cluster counts with newly assignment variant
            scs = scs.set_value(best_clus_array_idx, 'size', scs['size'].values[best_clus_pos]+1)

    return(scs, probs, ccert)

def collate_variant_output(var_df, var_filt_df, Nvar, ccert, probs, clus_probs, phis, scs, best_clus_list, snvs):
    lolim = 0 if snvs else 1
    uplim = 2 if snvs else 7
    nvar_filt = len(var_filt_df)
    var_df.index = range(nvar_filt, nvar_filt + Nvar)

    hpd_lo = ccert.columns.values[-2]
    hpd_hi = ccert.columns.values[-1]
    ccert_add = var_df[var_df.columns.values[lolim:uplim]]
    ccert_add['most_likely_assignment'] = best_clus_list
    ccert_add['average_ccf'] = np.array(phis)
    ccert_add[hpd_lo] = float('nan')
    ccert_add[hpd_hi] = float('nan')
    ccert = pd.concat([ccert, ccert_add], axis=0)
    ccert.index = range(len(ccert)) #reindex

    # NOTE: clus probabilities are of a different format here
    # CN likelihoods are converted to probabilities,
    # rather than probs being based on MCMC traces
    probs_add = var_df[var_df.columns.values[lolim:uplim]]
    uplim = uplim if snvs else uplim - 1

    prob_cols = probs.columns.values[uplim:]
    for idx,pc in enumerate(prob_cols):
        probs_add[pc] = 0.

    pa_clus_cols = ['cluster%d' % cid for cid in scs.clus_id.values]
    for idx,pac in enumerate(pa_clus_cols):
        probs_add[pac] = [round(cp, 2) for cp in clus_probs[:,idx]]

    probs = pd.concat([probs, probs_add], axis=0)
    probs.index = range(len(probs)) #reindex

    df = pd.concat([var_filt_df, var_df], axis=0)
    df.index = range(len(df)) #reindex

    return(df, ccert, probs)

def post_assign_vars(var_df, var_filt_df, rundir, sample, sparams, cparams, clus_th, snvs = False):
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

    # convert ccert to CCF values (expected input for write output function
    hpd_lo = 'HPD_lo'
    hpd_hi = 'HPD_hi'
    if len(ccert) > 0:
        hpd_lo = ccert.columns.values[-2]
        hpd_hi = ccert.columns.values[-1]
        ccert[hpd_lo] = ccert[hpd_lo]/purity
        ccert[hpd_hi] = ccert[hpd_hi]/purity

    filt_ids = get_var_ids(var_filt_df, snvs)
    ccert_ids = get_var_ids(ccert, snvs)

    if len(ccert) != len(var_filt_df) and snvs:
        # discrepancy due to subsampling
        var_df, var_filt_df = fix_subsampling_discrepancy(var_df, var_filt_df, filt_ids, ccert_ids)
        filt_ids = get_var_ids(var_filt_df, snvs)

    if len(ccert) != len(var_filt_df):
        raise ValueError('''Number of filtered variants does not match cluster output,
                         did you provide the correct filtered variants file?''')

    vartype = 'SNVs' if snvs else 'SVs'
    print('Post assigning %d %s...' % (len(var_df), vartype))

    # match ccert order
    matched_order = np.array([np.where(np.array(filt_ids)==x)[0][0] for x in ccert_ids])
    var_filt_df.index = range(len(var_filt_df))
    var_filt_df = var_filt_df.loc[matched_order]

    if snvs:
        if len(var_df) > 0:
            sup, dep, cn_states, Nvar, norm = load_data.get_snv_vals(var_df, male, cparams)
        if len(var_filt_df) > 0:
            fsup, fdep, fcn_states, fNvar, fnorm = load_data.get_snv_vals(var_filt_df, male, cparams)
    else:
        if len(var_df) > 0:
            sup, dep, cn_states, Nvar, norm = load_data.get_sv_vals(var_df, adjusted, male, cparams)
        if len(var_filt_df) > 0:
            fsup, fdep, fcn_states, fNvar, fnorm = load_data.get_sv_vals(var_filt_df, adjusted, male, cparams)

    if Nvar:
        no_cn_state = np.array([len(cn)==0 for cn in cn_states])
        if np.any(no_cn_state):
            keep = np.invert(no_cn_state)
            sup, dep, norm, cn_states = sup[keep], dep[keep], np.array(norm)[keep], cn_states[keep]
            var_df = var_df[keep]
            Nvar = Nvar - sum(no_cn_state)

    scs, orig_scs, best_clus_list, clus_probs, phis = filter_clusters(scs, clus_th, Nvar, \
                                                        sup, dep, norm, cn_states, purity)

    if len(scs) != len(orig_scs) and len(var_filt_df) > 0:
        scs, probs, ccert = reclassify_filtered_out_vars(scs, orig_scs, probs, ccert, \
                                                         fsup, fdep, fnorm, fcn_states, purity)

    # increment variant cluster counts with new variants assigned
    clusts = scs.clus_id.values
    for clus_idx in best_clus_list:
        clus = clusts[clus_idx]
        idx = scs.index[scs.clus_id.values==clus].values[0]
        scs = scs.set_value(idx, 'size', scs['size'][idx]+1)

    # convert from position index to cluster index (may differ)
    best_clus_list = [clusts[clus_idx] for clus_idx in best_clus_list]

    # collate variants from original output + post assigned
    if Nvar > 0:
        df, ccert, probs = collate_variant_output(var_df, var_filt_df, Nvar, ccert, probs, \
                                                  clus_probs, phis, scs, best_clus_list, snvs)
        sup = np.append(fsup, sup) if len(var_filt_df) > 0 else sup
        dep = np.append(fdep, dep) if len(var_filt_df) > 0 else dep
        norm = np.append(fnorm, norm) if len(var_filt_df) > 0 else norm
        cn_states = np.append(fcn_states, cn_states) if len(var_filt_df) > 0 else cn_states
    else:
        sup, dep, norm, cn_states = fsup, fdep, fnorm, fcn_states
        df = var_filt_df

    # retrieve traces for input to write_output func
    z_trace = np.loadtxt('%s/z_trace.txt.gz' % rundir)
    phi_trace = np.loadtxt('%s/phi_trace.txt.gz' % rundir)
    z_phi = run_clus.get_per_variant_phi(z_trace, phi_trace)

    if len(var_filt_df) < len(z_phi):
        # indicates coclustering, adjust z_phi length
        z_phi = z_phi[:len(var_filt_df)] if snvs else z_phi[len(z_phi)-len(var_filt_df):]

    z_phi_add = np.empty(len(var_df))
    z_phi_add[:] = np.NAN
    z_phi = np.concatenate([z_phi, z_phi_add])

    var_outfile = '%s/../%s_filtered_%s_post_assign.tsv' % (rundir, sample, vartype.lower())
    df.to_csv(var_outfile, sep='\t', index=False, na_rep='')

    print('Writing output to %s for %s' % (pa_outdir, vartype))
    clus_members = np.array([np.where(ccert.most_likely_assignment.values==i)[0] for i in scs.clus_id.values])
    write_output.write_out_files(df, scs, clus_members, probs, ccert,
                                 pa_outdir, sample, purity, sup, dep,
                                 norm, cn_states, fit, False, cnv_pval, z_phi, are_snvs = snvs)

def amend_coclus_results(rundir, sample, sparams):
    '''
    correct variant numbers in subclones output file by
    adding post assigned variants SVs and SNVs together
    '''
    purity   = sparams['pi']

    pa_outdir = '%s_post_assign/' % rundir
    sv_scs, sv_ccert, sv_probs, fit = load_data.get_run_output(pa_outdir, sample, purity, snvs = False)
    snv_scs, snv_ccert, snv_probs, fit = load_data.get_run_output(pa_outdir, sample, purity, snvs = True)

    scs = sv_scs.copy()
    for idx,sc in sv_scs.iterrows():
        sv_total = sum(sv_ccert.most_likely_assignment.values == int(sc.clus_id))
        snv_total = sum(snv_ccert.most_likely_assignment.values == int(sc.clus_id))
        scs = scs.set_value(idx, 'size', snv_total + sv_total)

    scs['CCF'] = scs.phi.values
    scs['phi'] = scs.phi.values * purity
    scs = scs[['clus_id','size','phi','CCF']]
    rename_cols =  {'clus_id': 'cluster', 'size': 'n_ssms', 'phi': 'proportion'}
    scs = scs.rename(columns = rename_cols)

    scs.to_csv('%s/%s_subclonal_structure.txt' % (pa_outdir, sample), sep='\t', index=False)
    scs.to_csv('%s/snvs/%s_subclonal_structure.txt' % (pa_outdir, sample), sep='\t', index=False)

def string_to_bool(v):
  return v.lower() in ("yes", "true", "t", "1")

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

    strict_cf = string_to_bool(Config.get('FilterParameters', 'strict_cnv_filt'))
    sv_offset = int(Config.get('FilterParameters', 'sv_offset'))
    gl_th     = int(Config.get('FilterParameters', 'germline_threshold'))
    cp_th     = float(Config.get('PostAssignParameters', 'clus_percent_threshold'))
    ca_th     = int(Config.get('PostAssignParameters', 'clus_absolute_threshold'))
    clus_th   = {'percent': cp_th, 'absolute': ca_th}

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

    if run.endswith("snvs"):
        sv_file=''

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

        if len(sv_df)>0:
            print('Matching copy-numbers for SVs...')
            sv_df = filt.match_copy_numbers(sv_df,cnv_df,strict_cf,sv_offset)
            sv_df = filt.match_copy_numbers(sv_df,cnv_df,strict_cf,sv_offset,\
                    ['chr2','pos2','dir2','classification','pos1'],'gtype2')

        if len(snv_df)>0:
            print('Matching copy-numbers for SNVs...')
            snv_df = filt.match_snv_copy_numbers(snv_df,cnv_df)
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

    copied_dir = False
    pa_outdir = '%s_post_assign/' % rundir
    if len(sv_df) > 0:
        sv_to_assign = get_var_to_assign(sv_df, sv_filt_df)
        sv_to_assign = filt.adjust_sv_read_counts(sv_to_assign, purity, ploidy, 0, rlen, Config)

        post_assign_vars(sv_to_assign, sv_filt_df, rundir, sample, sample_params,
                         cluster_params, clus_th)

    if len(snv_df) > 0:
        snv_to_assign = get_var_to_assign(snv_df, snv_filt_df, snvs = True)

        post_assign_vars(snv_to_assign, snv_filt_df, rundir, sample, sample_params,
                         cluster_params, clus_th, snvs = True)

    if len(sv_df) > 0 and len(snv_df) > 0:
        amend_coclus_results(rundir, sample, sample_params)
