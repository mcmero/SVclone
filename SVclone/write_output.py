import numpy as np
import pandas as pd
import os
import gzip
import itertools
import math
from operator import methodcaller

from . import dtypes
from . import cluster
from . import load_data

def adjust_vafs(mlcn_vect, ccert, vafs, pi, norm):
    '''
    adjust VAFs to CCFs using pv per variant
    '''
    mlcn  = pd.DataFrame(mlcn_vect)
    phis  = ccert.average_proportion.values / float(pi)
    probs = mlcn.pv.values / phis

    adj_vafs = (1 / probs) * vafs
    adj_vafs = np.array([av if av < 2 else 2 for av in adj_vafs])

    return(adj_vafs)

def dump_trace(center_trace, outf):
    outf = '%s.gz' % outf
    with gzip.open(outf, 'wt') as fout:
        df_traces = pd.DataFrame(center_trace)
        df_traces.to_csv(fout, sep='\t', index=False, header=False)

def write_out_files(df, clus_info, clus_members, df_probs, clus_cert, clus_out_dir, sample, pi, sup, dep, norm, cn_states, run_fit, smc_het, cnv_pval, z_phi, are_snvs=False):
    if are_snvs:
        clus_out_dir = '%s/snvs'%clus_out_dir
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)

    if len(run_fit) > 0:
        run_fit.to_csv('%s/%s_fit.txt' % (clus_out_dir, sample), sep='\t', index=False, header=False)

    clus_info['CCF'] = clus_info.phi.values
    clus_info['phi'] = clus_info.phi.values*pi #convert to proportions
    clus_info = clus_info[['clus_id','size','phi','CCF']]
    rename_cols =  {'clus_id': 'cluster', 'size': 'n_ssms', 'phi': 'proportion'}
    clus_info = clus_info.rename(columns = rename_cols)

    if smc_het:
        purity = min(max(pi, 0), 1)
        out_file = '%s/smc_1A_cellularity.txt' % clus_out_dir
        with open(out_file, 'w') as outf:
            outf.write('%f\n' % purity)
        with open('%s/smc_1B_number_of_clusters.txt' % clus_out_dir,'w') as outf:
            outf.write('%d\n' % len(clus_info))

        # convert cluster numbers to clean 1..n sequence
        smc_clus_info = clus_info.copy()
        smc_clus_info['new_clus'] = range(1, len(smc_clus_info)+1)

        mla = clus_cert[['most_likely_assignment']].copy()
        smc_cert = mla.copy()
        for idx, clus in smc_clus_info.iterrows():
            pd.options.mode.chained_assignment = None
            smc_cert.loc[mla.most_likely_assignment == clus.cluster, 'most_likely_assignment'] = int(clus.new_clus)
        smc_cert.to_csv('%s/smc_2A_mutations_to_clusters.txt' % clus_out_dir, sep='\t', index=False, header=False)

        smc_clus_info.cluster = smc_clus_info.new_clus
        smc_clus_info = smc_clus_info[['cluster', 'n_ssms', 'proportion']]
        smc_clus_info.to_csv('%s/smc_1C_cluster_structure.txt' % clus_out_dir, sep='\t', index=False, header=False)

    clus_info.to_csv('%s/%s_subclonal_structure.txt' % (clus_out_dir, sample), sep='\t', index=False)
    with open('%s/number_of_clusters.txt' % clus_out_dir,'w') as outf:
        outf.write("sample\tclusters\n")
        outf.write('%s\t%d\n'%(sample, len(clus_info)))

    cn_dtype    = dtypes.sv_cn_dtype    if not are_snvs else dtypes.snv_cn_dtype
    mlcn_dtype  = dtypes.sv_mlcn_dtype  if not are_snvs else dtypes.snv_mlcn_dtype
    mult_dtype  = dtypes.sv_mult_dtype  if not are_snvs else dtypes.snv_mult_dtype

    # re-sort all output in cluster member order
    cmem      = np.hstack(clus_members)
    df        = df.loc[cmem]
    df.index  = range(len(df)) #reindex df
    df_probs  = df_probs.loc[cmem]
    clus_cert = clus_cert.loc[cmem]
    z_phi     = z_phi[cmem]
    sup, dep, norm, cn_states = sup[cmem], dep[cmem], np.array(norm)[cmem], cn_states[cmem]

    # prepare data structures for output
    cn_vect     = np.empty((0, len(cmem)), dtype=cn_dtype)
    mlcn_vect   = np.empty((0, len(cmem)), dtype=mlcn_dtype)
    mult_vect   = np.empty((0, len(cmem)), dtype=mult_dtype)

    phis = clus_cert.average_ccf.values
    cns, pvs = cluster.get_most_likely_cn_states(cn_states, sup, dep, phis, pi, cnv_pval, norm)

    for idx, var in df.iterrows():
        gtype1, gtype2 = None, None
        if not are_snvs:
            if math.isnan(var.ID):
                continue
            gtype1, gtype2 = var['gtype1'].split('|'), var['gtype2'].split('|')
            gtype1, gtype2 = list(map(methodcaller('split', ','), gtype1)), list(map(methodcaller('split', ','), gtype2))

            #select the first clonal/major fraction copy-number state as the one to output
            maj_cn1, min_cn1 = [float(x) for x in gtype1[0]][:2] if gtype1[0][0]!='' else [0., 0.]
            maj_cn2, min_cn2 = [float(x) for x in gtype2[0]][:2] if gtype2[0][0]!='' else [0., 0.]

            tot_cn1 = maj_cn1 + min_cn1
            tot_cn2 = maj_cn2 + min_cn2 if not are_snvs else 0

            chr1 = str(var['chr1'])
            pos1 = int(var['pos1'])
            dir1 = str(var['dir1'])
            chr2 = str(var['chr2'])
            pos2 = int(var['pos2'])
            dir2 = str(var['dir2'])

            ref_cn, sc_cn, freq, frac = cns[idx]
            pv = pvs[idx]
            chrs_bearing_mut = int(sc_cn*freq) if (not np.isnan(sc_cn) or not np.isnan(freq)) else 0

            ref_cn = 0 if np.isnan(sc_cn) else ref_cn
            sc_cn = 0 if np.isnan(sc_cn) else sc_cn
            freq = 0 if np.isnan(freq) else freq

            cn_new_row = np.array([(chr1, pos1, dir1, tot_cn1, chrs_bearing_mut,
                                    chr2, pos2, dir2, tot_cn2, chrs_bearing_mut)], dtype=cn_dtype)
            pv_dev = np.abs(pv - var.adjusted_vaf)
            ml_new_row = np.array([(chr1, pos1, dir1, chr2, pos2, dir2, var['gtype1'], var['gtype2'],
                                    ref_cn, sc_cn, freq, sup[idx], dep[idx], frac, pv, pv_dev)],
                                    dtype=mlcn_dtype)

            var_states = cn_states[idx]
            tot_opts = ','.join(list(map(str,[int(x[1]) for x in var_states])))
            var_opts = ','.join(list(map(str,[int(round(x[1]*x[2])) for x in var_states])))
            probs =  cluster.get_probs(var_states, sup[idx], dep[idx], phis[idx], pi)
            
            mult_new_row = np.array([(chr1, pos1, dir1, chr2, pos2, dir2, sc_cn,
                                     int(round(freq*sc_cn)), tot_opts, var_opts, probs)], dtype=mult_dtype)

            cn_vect   = np.append(cn_vect, cn_new_row)
            mlcn_vect = np.append(mlcn_vect, ml_new_row)
            mult_vect = np.append(mult_vect, mult_new_row)
        else:
            gtype1 = list(map(lambda x: x.split(','), var['gtype'].split('|')))
            maj_cn1, min_cn1 = [float(x) for x in gtype1[0]][:2] if gtype1[0][0]!='' else [0., 0.]
            tot_cn1 = maj_cn1 + min_cn1

            chrom = str(var['chrom'])
            pos = int(var['pos'])

            ref_cn, sc_cn, freq, frac = cns[idx]
            pv = pvs[idx]
            pv_dev = np.abs(pv - float(var['var']) / (var['var'] + var.ref))

            cn_new_row = np.array([(chrom, pos, tot_cn1, int(sc_cn*freq))], dtype=cn_dtype)
            ml_new_row = np.array([(chrom, pos, var['gtype'], ref_cn, sc_cn, freq, frac, pv, pv_dev)],
                                  dtype=mlcn_dtype)

            var_states = cn_states[idx]
            tot_opts = ','.join(list(map(str,[int(x[1]) for x in var_states])))
            var_opts = ','.join(list(map(str,[int(round(x[1]*x[2])) for x in var_states])))
            probs =  cluster.get_probs(var_states, sup[idx], dep[idx], phis[idx], pi)

            mult_new_row = np.array([(chrom, pos, sc_cn, int(round(freq*sc_cn)), tot_opts, var_opts, probs)],
                    dtype=mult_dtype)

            cn_vect = np.append(cn_vect, cn_new_row)
            mlcn_vect = np.append(mlcn_vect, ml_new_row)
            mult_vect = np.append(mult_vect, mult_new_row)

    #adjust cluster freqs to cell prevalence
    hpd_lo = clus_cert.columns.values[-2]
    hpd_hi = clus_cert.columns.values[-1]
    clus_cert.average_ccf = clus_cert.average_ccf.values*pi
    clus_cert[hpd_lo] = clus_cert[hpd_lo].values*pi
    clus_cert[hpd_hi] = clus_cert[hpd_hi].values*pi

    #rename cols
    rename_cols =  {'average_ccf': 'average_proportion'}
    clus_cert = clus_cert.rename(columns = rename_cols)
    df_probs = df_probs.rename(columns = rename_cols)

    ccf_out = pd.DataFrame()
    if are_snvs:
        ccf_out = df[df.columns.values[:2]]
        vafs = df['var'].values / (df['var'].values + df['ref'].values)
        ccf_out['raw_VAF'] = vafs
        ccf_out['variant_CCF'] = adjust_vafs(mlcn_vect, clus_cert, vafs, pi, norm)
        ccf_out['variant_CCF_trace'] = z_phi
        ccf_out['cluster_proportion'] = clus_cert.average_proportion.values
        ccf_out['cluster_CCF'] = clus_cert.average_proportion.values / pi
    else:
        ccf_out = df[df.columns.values[:7]]
        ccf_out['raw_mean_VAF'] = df['raw_mean_vaf'].values
        ccf_out['adjusted_VAF'] = df['adjusted_vaf'].values
        ccf_out['variant_CCF'] = adjust_vafs(mlcn_vect, clus_cert, df['adjusted_vaf'].values, pi, norm)
        ccf_out['variant_CCF_trace'] = z_phi
        ccf_out['cluster_proportion'] = clus_cert.average_proportion.values
        ccf_out['cluster_CCF'] = clus_cert.average_proportion.values / pi

    pd.DataFrame(mlcn_vect).to_csv('%s/%s_most_likely_copynumbers.txt'%(clus_out_dir, sample), sep='\t', index=False)
    pd.DataFrame(mult_vect).to_csv('%s/%s_multiplicity.txt'%(clus_out_dir, sample), sep='\t', index=False)
    pd.DataFrame(cn_vect).to_csv('%s/%s_copynumber.txt'%(clus_out_dir, sample), sep='\t', index=False)
    df_probs.to_csv('%s/%s_assignment_probability_table.txt'%(clus_out_dir, sample), sep='\t', index=False)
    clus_cert.to_csv('%s/%s_cluster_certainty.txt'%(clus_out_dir, sample), sep='\t', index=False)
    ccf_out.to_csv('%s/%s_vaf_ccf.txt' % (clus_out_dir, sample), sep='\t', index=False)
