from __future__ import print_function

import pandas as pd
import numpy as np
import vcf
import configparser
import os
from SVclone.SVprocess import svp_load_data as svp_load

def get_normal_copynumber(chrom, male):
    if male and (chrom == 'X' or chrom == 'chrX'):
        return 1.
    elif chrom == 'Y' or chrom == 'chrY':
        return 1.
    else:
        return 2.

def get_sv_vals(sv_df, cparams):
    adjusted = cparams['adjusted']
    male = cparams['male']
    chroms = sv_df.chr1.values # consider only chrom 1 for normal copy-number
    norm_cn = [get_normal_copynumber(c, male) for c in chroms]
    norm1, norm2 = sv_df.adjusted_norm1.values, sv_df.adjusted_norm2.values
    norm1 = [int(round(n)) for n in norm1]
    norm2 = [int(round(n)) for n in norm2]

    Nvar = len(sv_df)
    sup = sv_df.adjusted_support.map(float).values if adjusted else sv_df.support.map(float).values

    return sup, norm1, norm2, Nvar, norm_cn

def get_snv_vals(df, cparams):
    male = cparams['male']
    n = df['ref'].map(float).values
    b = df['var'].map(float).values

    norm_cn = [get_normal_copynumber(c, male) for c in df.chrom.values]
    return b, (n+b), len(b), norm_cn

def load_svs(sv_file):
    dat = pd.read_csv(sv_file,delimiter='\t',dtype=None, low_memory=False)
    sv_df = pd.DataFrame(dat)
    sv_df.chr1 = sv_df.chr1.map(str)
    sv_df.chr2 = sv_df.chr2.map(str)
    sv_df['raw_norm_mean'] = [np.mean(norm) for norm in zip(sv_df['norm1'].values, sv_df['norm2'].values)]
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
        col_names = ['chr', 'startpos', 'endpos', 'norm_total', 'norm_minor', 'tumour_total', 'tumour_minor']
        cnv_df = pd.DataFrame(pd.read_csv(cnv_file, delimiter=',',
                              dtype=None, names=col_names, index_col=0,
                              skip_blank_lines=True))

    try:
        if 'nMaj1_A' in cnv_df.columns.values:
            # battenberg input
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.nMaj1_A.values))]
            cnv_df['chr'] = cnv_df['chr'].map(str)
            cnv_df['nMaj1_A'] = cnv_df['nMaj1_A'].map(float)
            cnv_df['nMin1_A'] = cnv_df['nMin1_A'].map(float)

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
            select_cols = ['chr', 'startpos', 'endpos', 'gtype']
            return cnv_df[select_cols]

        elif 'battenberg_nMaj1_A' in cnv_df.columns.values:
            # also battenberg
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.battenberg_nMaj1_A.values))]
            cnv_df['chr'] = cnv_df['chr'].map(str)
            cnv_df['battenberg_nMaj1_A'] = cnv_df['battenberg_nMaj1_A'].map(float)
            cnv_df['battenberg_nMin1_A'] = cnv_df['battenberg_nMin1_A'].map(float)

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
            select_cols = ['chr', 'startpos', 'endpos', 'gtype']
            return cnv_df[select_cols]

        elif 'clonal_frequency' in cnv_df.columns.values:
            # pcawg clonal copy-number calls format
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.major_cn.values))]
            cnv_df['chromosome'] = cnv_df.chromosome.map(str)
            cnv_df = cnv_df.rename(columns={'chromosome': 'chr', 'start': 'startpos', 'end': 'endpos'})
            gtypes = cnv_df.major_cn.map(str) + ',' \
                        + cnv_df.minor_cn.map(str) \
                        + ',' + cnv_df.clonal_frequency.map(str)
            cnv_df['gtype'] = gtypes
            select_cols = ['chr', 'startpos', 'endpos', 'gtype']
            return cnv_df[select_cols]

        elif 'star' in cnv_df.columns.values:
            # pcawg star copy-number calls format
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.major_cn.values))]
            cnv_df['chromosome'] = cnv_df.chromosome.map(str)
            cnv_df = cnv_df.rename(columns={'chromosome': 'chr', 'start': 'startpos', 'end': 'endpos'})
            gtypes = cnv_df.major_cn.map(str) + ',' + cnv_df.minor_cn.map(str) + ',1.0'
            cnv_df['gtype'] = gtypes
            select_cols = ['chr', 'startpos', 'endpos', 'gtype']
            return cnv_df[select_cols]

        else:
            # caveman input
            cnv_df = cnv_df[np.invert(np.isnan(cnv_df.tumour_total.values))]
            major = cnv_df.tumour_total.map(int) - cnv_df.tumour_minor.map(int)
            gtypes = major.map(str) + ',' + cnv_df.tumour_minor.map(str) + ',1.0'
            cnv_df['gtype'] = gtypes
            select_cols = ['chr', 'startpos', 'endpos', 'gtype']
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
    snv_out = snv_out[['chrom', 'pos', 'gtype', 'ref', 'var']]
    return snv_out

def load_snvs_multisnv(snvs, sample):
    snv_dtype = [('chrom', '<U50'), ('pos', int), ('gtype', '<U50'), ('ref', float), ('var', float)]
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
    snv_dtype = [('chrom', '<U50'), ('pos',int), ('gtype', '<U50'), ('ref',float), ('var',float)]
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
    snv_dtype = [('chrom', '<U50'), ('pos',int), ('gtype', '<U50'), ('ref',float), ('var',float)]
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
            tmp = np.array((record.CHROM,record.POS, '',ref_reads,variant_reads),dtype=snv_dtype)
            snv_df = np.append(snv_df,tmp)

    return pd.DataFrame(snv_df)

def load_snvs_sanger(snvs):
    vcf_reader = vcf.Reader(filename=snvs)
    snv_dtype = [('chrom', '<U50'), ('pos',int), ('gtype', '<U50'), ('ref',float), ('var',float)]
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
            tmp = np.array((record.CHROM,record.POS, '', ref_reads,variant_reads),dtype=snv_dtype)
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

    Config = configparser.ConfigParser()
    cfg_file = Config.read(cfg)

    if len(cfg_file)==0:
        raise ValueError('No configuration file found')

    mean_cov        = float(Config.get('BamParameters', 'mean_cov'))
    clus_limit      = int(Config.get('ClusterParameters', 'clus_limit'))
    n_iter          = int(Config.get('ClusterParameters', 'n_iter'))
    threads         = int(Config.get('ClusterParameters', 'threads'))
    repeat          = int(Config.get('ClusterParameters', 'repeat'))
    adjusted        = string_to_bool(Config.get('ClusterParameters', 'adjusted'))
    male            = string_to_bool(Config.get('ClusterParameters', 'male'))
    sv_to_sim       = int(Config.get('ClusterParameters', 'sv_to_sim'))

    if XX:
        male = False
    if XY:
        male = True

    pi, pl = svp_load.get_purity_ploidy(pp_file, sample, out)
    rlen, insert, std = svp_load.get_read_params(param_file, sample, out)

    sample_params  = { 'sample': sample, 'ploidy': pl, 'pi': pi, 'rlen': rlen,
                       'insert': insert, 'mean_cov': mean_cov }
    cluster_params = { 'n_iter': n_iter, 'male': male, 'adjusted': adjusted,
                       'sv_to_sim': sv_to_sim, 'threads': threads, 'repeat': repeat,
                       'clus_limit': clus_limit }

    return sample_params, cluster_params
