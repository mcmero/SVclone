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

def get_ll_probs(s, d, n, cn, purity, scs):
    lls = np.array(np.max(cluster.calc_lik(cn, s, d, scs.phi.values[0], purity, n)[1]))
    for j in range(1, len(scs)):
        ll = cluster.calc_lik(cn, s, d, scs.phi.values[j], purity, n)[1]
        lls = np.append(lls, np.max(ll))
    try:
        lls_conv = [math.e ** ll for ll in lls]
        ll_probs = [float(lc) / sum(lls_conv) for lc in lls_conv]
    except TypeError:
        ll_probs = (math.e ** lls) / lls

    return lls, ll_probs

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

    if snvs:
        sup, dep, cn_states, Nvar, norm = load_data.get_snv_vals(var_df, male)
        if len(var_filt_df) > 0:
            fsup, fdep, fcn_states, fNvar, fnorm = load_data.get_snv_vals(var_filt_df, male)
    else:
        sup, dep, cn_states, Nvar, norm = load_data.get_sv_vals(var_df, adjusted, male)
        if len(var_filt_df) > 0:
            fsup, fdep, fcn_states, fNvar, fnorm = load_data.get_sv_vals(var_filt_df, adjusted, male)

    if len(ccert) != len(var_filt_df):
        raise ValueError('''Number of filtered variants does not match cluster output,
                         did you provide the correct filtered variants file?''')

    no_cn_state = np.array([len(cn)==0 for cn in cn_states])
    if np.any(no_cn_state):
        keep = np.invert(no_cn_state)
        sup, dep, norm, cn_states = sup[keep], dep[keep], np.array(norm)[keep], cn_states[keep]
        var_df = var_df[keep]
        Nvar = Nvar - sum(no_cn_state)

    best_clus_list = np.array([])
    clus_probs = np.array([])

    # filter out low proportion clusters
    percent_filt = (scs['size']/sum(scs['size'])).values > clus_th['percent']
    abs_filt     = scs['size'].values > clus_th['absolute']
    orig_scs     = scs.copy()
    scs          = scs[np.logical_and(percent_filt, abs_filt)]

    if len(scs) > 1:
        for i in range(Nvar):
            lls, ll_probs = get_ll_probs(sup[i], dep[i], norm[i], cn_states[i], purity, scs)
            clus_probs = np.concatenate([clus_probs, [ll_probs]], axis=0) if len(clus_probs)>0 else np.array([ll_probs])
            best_clus = np.where(np.max(lls)==lls)[0][0]
            best_clus_list = np.append(best_clus_list, [best_clus])
    else:
        best_clus_list = [0] * Nvar
        clus_probs     = np.array([[1.]] * Nvar)

    best_clus_list = np.array(map(int, best_clus_list))
    phis = scs.phi.values[best_clus_list]

    # reclassify minor cluster variants if filtered out
    if len(scs) != len(orig_scs) and len(var_filt_df) > 0:
        scs_fo = orig_scs[np.invert(np.logical_and(percent_filt, abs_filt))]
        # drops cols from probabilities dataframe
        for i in scs_fo.clus_id.values:
            probs = probs.drop('cluster%d' % i, 1)

        mla = ccert.most_likely_assignment.values
        which_reclass = np.where(np.array([cv in scs_fo.clus_id.values for cv in mla]))[0]
        if len(which_reclass) > 0:
            print('Reassigning %d variants from filtered out clusters' % len(which_reclass))
            for i in which_reclass:
                # find most likely cluster for variant
                lls, ll_probs = get_ll_probs(fsup[i], fdep[i], fnorm[i], fcn_states[i], purity, scs)
                best_clus = np.where(np.max(lls)==lls)[0][0]
                best_clus_id = scs.index[scs.clus_id.values==scs.clus_id.values[best_clus]].values[0]

                # change variant's most likely assignment to new membership
                ccert = ccert.set_value(i, 'most_likely_assignment', best_clus_id)

                # increment subclonal cluster counts with newly assignment variant
                idx = scs.index[scs.clus_id.values==best_clus_id].values[0]
                scs = scs.set_value(idx, 'size', scs['size'][idx]+1)

    clusts = scs.clus_id.values
    for clus_idx in best_clus_list:
        clus = clusts[clus_idx]
        idx = scs.index[scs.clus_id.values==clus].values[0]
        scs = scs.set_value(idx, 'size', scs['size'][idx]+1)

    # convert from position index to cluster index (may differ)
    best_clus_list = [clusts[clus_idx] for clus_idx in best_clus_list]

    hpd_lo = 'HPD_lo'
    hpd_hi = 'HPD_hi'
    if len(ccert) > 0:
        hpd_lo = ccert.columns.values[-2]
        hpd_hi = ccert.columns.values[-1]
        ccert[hpd_lo] = ccert[hpd_lo]/purity
        ccert[hpd_hi] = ccert[hpd_hi]/purity

    lolim = 0 if snvs else 1
    uplim = 2 if snvs else 7
    nvar, nvar_filt = len(var_df), len(var_filt_df)
    var_df.index = range(nvar_filt, nvar_filt + nvar)
    ccert_add = var_df[var_df.columns.values[lolim:uplim]]
    ccert_add['most_likely_assignment'] = best_clus_list
    ccert_add['average_ccf'] = np.array(phis)
    ccert_add[hpd_lo] = float('nan')
    ccert_add[hpd_hi] = float('nan')
    ccert = pd.concat([ccert, ccert_add], axis=0)

    z_trace = np.loadtxt('%s/z_trace.txt.gz' % rundir)
    phi_trace = np.loadtxt('%s/phi_trace.txt.gz' % rundir)
    z_phi = run_clus.get_per_variant_phi(z_trace, phi_trace)
    if len(var_filt_df) < len(z_phi):
        # indicates coclustering, adjust z_phi length
        z_phi = z_phi[:len(var_filt_df)] if snvs else z_phi[len(z_phi)-len(var_filt_df):]

    z_phi_add = np.empty(len(var_df))
    z_phi_add[:] = np.NAN
    z_phi = np.concatenate([z_phi, z_phi_add])

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

    copied_dir = False
    pa_outdir = '%s_post_assign/' % rundir
    if len(sv_df) > 0:
        sv_to_assign = get_var_to_assign(sv_df, sv_filt_df)
        if len(sv_to_assign) > 0:
            sv_to_assign = filt.adjust_sv_read_counts(sv_to_assign, purity, ploidy, 0, rlen, Config)
            print('Post assigning %d SVs...' % len(sv_to_assign))
            post_assign_vars(sv_to_assign, sv_filt_df, rundir, sample, sample_params,
                             cluster_params, clus_th)
        else:
            # copy best run dir (no need to reassign)
            # remove dir if it exists
            if os.path.exists(pa_outdir):
                rmtree(pa_outdir)
            copytree(rundir, pa_outdir, ignore=ignore_patterns('*.gz'))
            copied_dir = True

    if len(snv_df) > 0:
        snv_to_assign = get_var_to_assign(snv_df, snv_filt_df, snvs = True)
        if len(snv_to_assign) > 0:
            print('Post assigning %d SNVs...' % len(snv_to_assign))
            post_assign_vars(snv_to_assign, snv_filt_df, rundir, sample, sample_params,
                             cluster_params, clus_th, snvs = True)
        elif not copied_dir:
            # copy snvs dir (no need to post-assign)
            # remove dir if it exists
            if os.path.exists(pa_outdir):
                rmtree(pa_outdir)
	    if run == 'best_run_snvs':
            	copytree('%s/' % rundir, '%s_post_assign/snvs' % rundir, ignore=ignore_patterns('*.gz'))
            else:
		copytree('%s/snvs' % rundir, '%s_post_assign/snvs' % rundir, ignore=ignore_patterns('*.gz'))

    if len(sv_df) > 0 and len(snv_df) > 0:
        amend_coclus_results(rundir, sample, sample_params)
