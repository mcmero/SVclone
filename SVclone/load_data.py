from __future__ import print_function

import pandas as pd
import numpy as np
import vcf
import ConfigParser
import os

from . import cluster
from SVprocess import svp_load_data as svp_load

def get_normal_copynumber(chrom, male):
    if male and (chrom == 'X' or chrom == 'chrX'):
        return 1.
    elif chrom == 'Y' or chrom == 'chrY':
        return 1.
    else:
        return 2.

def get_sv_vals(sv_df, adjusted, male, cparams):
    combos = sv_df.apply(cluster.get_sv_allele_combos, axis=1, args=(cparams,)).values
    sides = sv_df.preferred_side.map(int).values
    cn_states = [cn[side] for cn,side in zip(combos, sides)]
    cn_states = pd.DataFrame([[sv] for sv in cn_states])[0].values
    chroms = [c1 if side == 0 else c2 for c1,c2,side in zip(sv_df.chr1.values, sv_df.chr2.values, sides)]
    norm = [get_normal_copynumber(c, male) for c in chroms]

    if adjusted:
        sup = sv_df.adjusted_support.map(float).values
        dep = sv_df.adjusted_depth.map(float).values
        Nvar = len(sv_df)
        return sup,dep,cn_states,Nvar,norm
    else:
        sup  = sv_df.support.map(float).values
        norm = zip(sv_df.norm1.values,sv_df.norm2.values)
        norm = np.array([float(n[side]) for n,side in zip(norm,sv_df.preferred_side.values)])
        dep  = norm+sup
        Nvar = len(sv_df)
        return sup,dep,cn_states,Nvar,norm

def get_snv_vals(df, male, cparams):
    n = df['ref'].map(float).values
    b = df['var'].map(float).values

    cn_states = [cluster.get_allele_combos(gtype.split('|'), cparams) for gtype in df.gtype]
    cn_states = pd.DataFrame([[cn] for cn in cn_states])[0].values

    norm = [get_normal_copynumber(c, male) for c in df.chrom.values]
    return b,(n+b),cn_states,len(b),norm

def load_svs(sv_file):
    dat = pd.read_csv(sv_file,delimiter='\t',dtype=None, low_memory=False)
    sv_df = pd.DataFrame(dat)
    sv_df.chr1 = sv_df.chr1.map(str)
    sv_df.chr2 = sv_df.chr2.map(str)
    sv_df['norm_mean'] = map(np.mean,zip(sv_df['norm1'].values,sv_df['norm2'].values))
    return sv_df

def load_cnvs(cnv_file):
    cnv = pd.read_csv(cnv_file,delimiter='\t',dtype=None)
    cnv_df = pd.DataFrame(cnv)
    if len(cnv_df) == 0:
        print('WARNING: Copy-number file contains no records.')
        cnv_df['chr'] = ''
        cnv_df['startpos'] = float('nan')
        cnv_df['endpos'] = float('nan')
        cnv_df['gtype'] = ''
        return cnv_df

    if len(cnv_df.columns)==1:
        # assume caveman csv file
        col_names = ['chr','startpos','endpos','norm_total','norm_minor','tumour_total','tumour_minor']
        cnv_df = pd.DataFrame(pd.read_csv(cnv_file,delimiter=',',dtype=None,names=col_names,index_col=0,skip_blank_lines=True))

    try:
        if 'nMaj1_A' in cnv_df.columns.values:
            # battenberg input
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.nMaj1_A.values))]
            cnv_df['chr'] = map(str,cnv_df['chr'])
            cnv_df['nMaj1_A'] = map(float,cnv_df['nMaj1_A'])
            cnv_df['nMin1_A'] = map(float,cnv_df['nMin1_A'])

            gtypes = cnv_df['nMaj1_A'].map(str) + ',' + \
                     cnv_df['nMin1_A'].map(str) + ',' + \
                     cnv_df['frac1_A'].map(str)

            # join subclonal genotypes
            subclonal = cnv_df['frac1_A']!=1
            cnv_sc    = cnv_df[subclonal]
            gtypes[subclonal] = gtypes[subclonal] + '|' + \
                                cnv_sc['nMaj2_A'].map(str) + ',' + \
                                cnv_sc['nMin2_A'].map(str) + ',' + \
                                cnv_sc['frac2_A'].map(str)

            cnv_df['gtype'] = gtypes
            select_cols = ['chr','startpos','endpos','gtype']
            return cnv_df[select_cols]

        elif 'battenberg_nMaj1_A' in cnv_df.columns.values:
            # also battenberg
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.battenberg_nMaj1_A.values))]
            cnv_df['chr'] = map(str,cnv_df['chr'])
            cnv_df['battenberg_nMaj1_A'] = map(float,cnv_df['battenberg_nMaj1_A'])
            cnv_df['battenberg_nMin1_A'] = map(float,cnv_df['battenberg_nMin1_A'])

            gtypes = cnv_df['battenberg_nMaj1_A'].map(str) + ',' + \
                     cnv_df['battenberg_nMin1_A'].map(str) + ',' + \
                     cnv_df['battenberg_frac1_A'].map(str)

            # join subclonal genotypes
            subclonal = cnv_df['battenberg_frac1_A']!=1
            cnv_sc    = cnv_df[subclonal]
            gtypes[subclonal] = gtypes[subclonal] + '|' + \
                                cnv_sc['battenberg_nMaj2_A'].map(str) + ',' + \
                                cnv_sc['battenberg_nMin2_A'].map(str) + ',' + \
                                cnv_sc['battenberg_frac2_A'].map(str)

            cnv_df['gtype'] = gtypes
            cnv_df = cnv_df.rename(columns={'start': 'startpos', 'end': 'endpos'})
            select_cols = ['chr','startpos','endpos','gtype']
            return cnv_df[select_cols]

        elif 'clonal_frequency' in cnv_df.columns.values:
            # pcawg clonal copy-number calls format
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.major_cn.values))]
            cnv_df['chromosome'] = cnv_df.chromosome.map(str)
            cnv_df = cnv_df.rename(columns={'chromosome': 'chr', 'start': 'startpos', 'end': 'endpos'})
            gtypes = cnv_df.major_cn.map(str) + ',' + cnv_df.minor_cn.map(str) + ',' + cnv_df.clonal_frequency.map(str)
            cnv_df['gtype'] = gtypes
            select_cols = ['chr','startpos','endpos','gtype']
            return cnv_df[select_cols]

        elif 'star' in cnv_df.columns.values:
            # pcawg star copy-number calls format
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.major_cn.values))]
            cnv_df['chromosome'] = cnv_df.chromosome.map(str)
            cnv_df = cnv_df.rename(columns={'chromosome': 'chr', 'start': 'startpos', 'end': 'endpos'})
            gtypes = cnv_df.major_cn.map(str) + ',' + cnv_df.minor_cn.map(str) + ',1.0'
            cnv_df['gtype'] = gtypes
            select_cols = ['chr','startpos','endpos','gtype']
            return cnv_df[select_cols]

        else:
            # caveman input
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.tumour_total.values))]
            major = cnv_df.tumour_total.map(int) - cnv_df.tumour_minor.map(int)
            gtypes = major.map(str) + ',' + cnv_df.tumour_minor.map(str) + ',1.0'
            cnv_df['gtype'] = gtypes
            select_cols = ['chr','startpos','endpos','gtype']
            return cnv_df[select_cols]

    except KeyError:
        raise Exception('''CNV file column names not recognised.
        Check the input is a battenberg output file (with headers) or an ASCAT caveman CSV.''')

def load_snvs_mutect_callstats(snvs):
    snv_df = pd.DataFrame(pd.read_csv(snvs,delimiter='\t',low_memory=False,dtype=None,comment="#"))
    snv_df = snv_df[snv_df['judgement']=='KEEP']

    snv_out = {'chrom' : snv_df.contig.map(str),
               'pos'   : snv_df.position.map(int),
               'gtype' : '',
               'ref'   : snv_df.t_ref_sum.map(float),
               'var'   : snv_df.t_alt_sum.map(float)}

    snv_out = pd.DataFrame(snv_out)
    snv_out = snv_out[['chrom','pos','gtype','ref','var']]
    return snv_out

def load_snvs_multisnv(snvs, sample):
    snv_dtype = [('chrom','S50'),('pos',int),('gtype','S50'),('ref',float),('var',float)]
    snv_df = np.empty([0,5],dtype=snv_dtype)

    vcf_reader = vcf.Reader(filename=snvs)
    samples = vcf_reader.samples
    if sample not in samples:
        raise Exception('VCF SNV file is of invalid format: sample %s not present in VCF samples.' % sample)

    for record in vcf_reader:
        if record.FILTER is not None:
            if len(record.FILTER)>0:
                continue

        allele_counts = record.genotype(sample)['BCOUNT']
        tumor_reads = {
            'A': allele_counts[0],
            'C': allele_counts[1],
            'G': allele_counts[2],
            'T': allele_counts[3],
        }

        ref_reads     = tumor_reads[record.REF]
        variant_reads = tumor_reads[str(record.ALT[0])]
        total_reads   = ref_reads + variant_reads

        if variant_reads != 0:
            tmp = np.array((record.CHROM, record.POS, '', ref_reads, variant_reads), dtype=snv_dtype)
            snv_df = np.append(snv_df, tmp)

    return pd.DataFrame(snv_df)

def load_snvs_consensus(snvs):
    vcf_reader = vcf.Reader(filename=snvs)
    snv_dtype = [('chrom','S50'),('pos',int),('gtype','S50'),('ref',float),('var',float)]
    snv_df = np.empty([0,5],dtype=snv_dtype)

    for record in vcf_reader:
        try:
            ref_reads, variant_reads = record.INFO['t_ref_count'], record.INFO['t_alt_count']
            total_reads = ref_reads + variant_reads
            if variant_reads != 0:
                tmp = np.array((record.CHROM, record.POS, '', ref_reads, variant_reads), dtype=snv_dtype)
                snv_df = np.append(snv_df,tmp)
        except KeyError:
            print('WARNING: missing count field(s) in record %s:%d' % (record.CHROM, record.POS))

    return pd.DataFrame(snv_df)

def load_snvs_mutect(snvs,sample):
    vcf_reader = vcf.Reader(filename=snvs)
    snv_dtype = [('chrom','S50'),('pos',int),('gtype','S50'),('ref',float),('var',float)]
    snv_df = np.empty([0,5],dtype=snv_dtype)

    samples = vcf_reader.samples
    if len(samples)==0:
        raise Exception('No samples found in VCF!')
    elif not np.any(np.array(samples)==sample):
        print('Warning, sample not found in VCF, selecting first sample')
        sample = samples[0]

    for record in vcf_reader:
        if record.FILTER is not None:
            if len(record.FILTER)>0:
                continue
        try:
            if record.genotype('normal')['AD'][1]>0:
                print('Removing variant %s:%d as it contains reads in the germline.' % (record.CHROM, record.POS))
                continue
        except KeyError:
            pass
        ad = record.genotype(sample)['AD']
        ref_reads, variant_reads = float(ad[0]), float(ad[1])
        total_reads = ref_reads + variant_reads
        if variant_reads!=0:
            tmp = np.array((record.CHROM,record.POS,'',ref_reads,variant_reads),dtype=snv_dtype)
            snv_df = np.append(snv_df,tmp)

    return pd.DataFrame(snv_df)

def load_snvs_sanger(snvs):
    vcf_reader = vcf.Reader(filename=snvs)
    snv_dtype = [('chrom','S50'),('pos',int),('gtype','S50'),('ref',float),('var',float)]
    snv_df = np.empty([0,5],dtype=snv_dtype)

    #code adapted from: https://github.com/morrislab/phylowgs/blob/master/parser/create_phylowgs_inputs.py
    samples = vcf_reader.samples
    if samples[0].lower()!='normal' or (samples[1].lower()!='tumour' and samples[1].lower()!='tumor'):
        raise Exception('VCF SNV file is of invalid format. Expected "NORMAL" and "TUMOUR" samples.')

    for record in vcf_reader:
        # get most likely genotypes
        genotypes = []
        broad_syn = False
        if 'TG' in record.INFO and 'SG' in record.INFO:
          genotypes = [record.INFO['TG'], record.INFO['SG']]
        else:
          genotypes = [x for x in record.INFO.keys() if "/" in x ]
          broad_syn = True

        if len(genotypes)==0:
            print('Warning: no valid genotypes for variant %s:%d; skipping.'%(record.CHROM,record.POS))
            continue
        if record.FILTER is not None:
            if len(record.FILTER)>0:
                continue

        variant_set = set()
        reference_nt = ''
        while len(variant_set) == 0:
            if len(genotypes) == 0:
                break
                #raise Exception('No more genotypes to find variant_nt in for %s' % variant)
            if broad_syn:
                gt = [x for x in genotypes if x.split("/")[0] != x.split("/")[1]][0]
                tumour_gt, normal_gt = gt.split('/')
            else:
                gt = genotypes.pop(0)
                normal_gt, tumour_gt = gt.split('/')
            if normal_gt[0] == normal_gt[1]:
                reference_nt = normal_gt[0]
                variant_set = set(tumour_gt) - set(reference_nt)
        variant_nt = variant_set.pop() if len(variant_set)!=0 else ''

        if variant_nt=='':
            print('Warning: no valid genotypes for variant %s:%d; skipping.'%(record.CHROM,record.POS))
            continue

        normal = record.genotype('NORMAL')
        tumour = record.genotype('TUMOUR')

        tumor_reads = {
          'forward': {
            'A': int(tumour['FAZ']),
            'C': int(tumour['FCZ']),
            'G': int(tumour['FGZ']),
            'T': int(tumour['FTZ']),
          },
          'reverse': {
            'A': int(tumour['RAZ']),
            'C': int(tumour['RCZ']),
            'G': int(tumour['RGZ']),
            'T': int(tumour['RTZ']),
          },
        }

        ref_reads = tumor_reads['forward'][reference_nt] + tumor_reads['reverse'][reference_nt]
        variant_reads = tumor_reads['forward'][variant_nt] + tumor_reads['reverse'][variant_nt]
        total_reads = ref_reads + variant_reads

        if variant_reads!=0:
            tmp = np.array((record.CHROM,record.POS,'',ref_reads,variant_reads),dtype=snv_dtype)
            snv_df = np.append(snv_df,tmp)

    #import matplotlib.pyplot as plt
    #fig, axes = plt.subplots(1, 1, sharex=False, sharey=False)
    #dep = snv_df['ref']+snv_df['var']
    #sup = map(float,snv_df['var'])
    #axes.hist(sup/dep);plt.savefig('/home/mcmero/Desktop/test')

    return pd.DataFrame(snv_df)

def string_to_bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def get_params_cluster_step(sample, cfg, out, pp_file, param_file, XX, XY):
    '''
    Load in paramaters used in cluster step from config file
    '''

    Config = ConfigParser.ConfigParser()
    cfg_file = Config.read(cfg)

    if len(cfg_file)==0:
        raise ValueError('No configuration file found')

    mean_cov        = float(Config.get('BamParameters', 'mean_cov'))
    shape           = float(Config.get('BetaParameters', 'alpha'))
    scale           = float(Config.get('BetaParameters', 'beta'))
    fixed_alpha     = Config.get('BetaParameters', 'fixed_alpha')
    phi_limit       = float(Config.get('ClusterParameters', 'phi_limit'))
    clus_limit      = int(Config.get('ClusterParameters', 'clus_limit'))
    subclone_diff   = float(Config.get('ClusterParameters', 'subclone_diff'))
    hpd_alpha       = float(Config.get('ClusterParameters', 'hpd_alpha'))

    n_runs          = int(Config.get('ClusterParameters', 'n_runs'))
    n_iter          = int(Config.get('ClusterParameters', 'n_iter'))
    burn            = int(Config.get('ClusterParameters', 'burn'))
    thin            = int(Config.get('ClusterParameters', 'thin'))
    threads         = int(Config.get('ClusterParameters', 'threads'))
    nclus_init      = Config.get('ClusterParameters', 'nclus_init')
    restrict_cnss   = string_to_bool(Config.get('ClusterParameters', 'restrict_cnv_search_space'))

    use_map         = string_to_bool(Config.get('ClusterParameters', 'map'))
    merge_clusts    = string_to_bool(Config.get('ClusterParameters', 'merge'))
    cocluster       = string_to_bool(Config.get('ClusterParameters', 'cocluster'))
    adjusted        = string_to_bool(Config.get('ClusterParameters', 'adjusted'))
    cnv_pval        = float(Config.get('ClusterParameters', 'clonal_cnv_pval'))
    adjust_phis     = string_to_bool(Config.get('ClusterParameters', 'adjust_phis'))
    male            = string_to_bool(Config.get('ClusterParameters', 'male'))
    sv_to_sim       = int(Config.get('ClusterParameters', 'sv_to_sim'))

    plot            = string_to_bool(Config.get('OutputParameters', 'plot'))
    ccf_reject      = float(Config.get('OutputParameters', 'ccf_reject_threshold'))
    smc_het         = string_to_bool(Config.get('OutputParameters', 'smc_het'))
    fit_metric      = Config.get('OutputParameters', 'fit_metric')
    cluster_penalty = int(Config.get('OutputParameters', 'cluster_penalty'))

    try:
        merge_iter      = int(Config.get('ClusterParameters', 'merge_iter'))
        merge_burn      = int(Config.get('ClusterParameters', 'merge_burn'))
    except ConfigParser.NoOptionError:
        merge_iter = int(round(n_iter / 4))
        merge_burn = int(round(burn / 4))

    if burn == 0 and use_map:
        print('No burn-in period specified, setting MAP to false.')
        use_map = False

    if XX:
        male = False
    if XY:
        male = True

    pi, pl = svp_load.get_purity_ploidy(pp_file, sample, out)
    rlen, insert, std = svp_load.get_read_params(param_file, sample, out)

    sample_params  = { 'sample': sample, 'ploidy': pl, 'pi': pi, 'rlen': rlen,
                       'insert': insert, 'mean_cov': mean_cov }
    cluster_params = { 'n_runs': n_runs, 'n_iter': n_iter, 'burn': burn, 'thin': thin, 'alpha': shape,
                       'beta': scale, 'use_map': use_map, 'hpd_alpha': hpd_alpha, 'fixed_alpha': fixed_alpha,
                       'male': male, 'merge_clusts': merge_clusts, 'adjusted': adjusted, 'phi_limit': phi_limit,
                       'clus_limit': clus_limit, 'subclone_diff': subclone_diff, 'cocluster': cocluster ,
                       'clonal_cnv_pval': cnv_pval, 'adjust_phis': adjust_phis, 'sv_to_sim': sv_to_sim,
                       'threads': threads, 'ccf_reject': ccf_reject, 'nclus_init': nclus_init,
                       'restrict_cnss': restrict_cnss, 'merge_iter': merge_iter, 'merge_burn': merge_burn }
    output_params  = { 'plot': plot, 'smc_het': smc_het, 'cluster_penalty': cluster_penalty, 'fit_metric': fit_metric }

    return sample_params, cluster_params, output_params

def get_run_output(rundir, sample, purity, snvs=False):
    outdir = rundir if not snvs or os.path.basename(rundir) == 'best_run_snvs' else '%s/snvs/' % rundir

    fit_file   = '%s/%s_fit.txt' % (outdir, sample)
    scs_file   = '%s/%s_subclonal_structure.txt' % (outdir, sample)
    ccert_file = '%s/%s_cluster_certainty.txt' % (outdir, sample)
    probs_file = '%s/%s_assignment_probability_table.txt' % (outdir, sample)

    fit   = pd.DataFrame()
    scs   = pd.DataFrame()
    ccert = pd.DataFrame()
    probs = pd.DataFrame()

    alt_scs_file = '%s/snvs/%s_subclonal_structure.txt' % (rundir, sample)
    alt_scs_file = alt_scs_file if not snvs else '%s/%s_subclonal_structure.txt' % (rundir, sample)

    if os.path.exists(scs_file):
        scs = pd.read_csv(scs_file, delimiter = '\t', dtype = None, low_memory = False)
    elif os.path.exists(alt_scs_file):
        scs = pd.read_csv(alt_scs_file, delimiter = '\t', dtype = None, low_memory = False)
    else:
        raise ValueError('No subclonal structure file exists!')

    # have to rename columns for function compatibility
    rename_cols =  {'cluster': 'clus_id', 'n_variants': 'size', 'CCF': 'phi'}
    if 'n_ssms' in scs.columns.values:
        rename_cols =  {'cluster': 'clus_id', 'n_ssms': 'size', 'CCF': 'phi'}
    scs = scs.rename(columns = rename_cols)
    scs = scs.drop('proportion', 1)

    if os.path.exists(ccert_file):
        ccert = pd.read_csv(ccert_file, delimiter = '\t', dtype = None, low_memory = False)

    if os.path.exists(probs_file):
        probs = pd.read_csv(probs_file, delimiter = '\t', dtype = None, low_memory = False)

    if os.path.exists(fit_file):
        fit = pd.read_csv(fit_file, delimiter = '\t', dtype = None, low_memory = False)

    if len(ccert) > 0:
        # have to rename columns for function compatibility
        rename_cols =  {'average_proportion': 'average_ccf'}
        ccert = ccert.rename(columns = rename_cols)
        ccert['average_ccf'] = ccert.average_ccf.values / purity

    return scs, ccert, probs, fit
