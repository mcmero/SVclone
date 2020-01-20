import nose2
import unittest
import numpy as np
import configparser
import os
import pandas as pd
import subprocess
from unittest import TestCase
from SVclone.SVprocess import bamtools
from SVclone.SVprocess import annotate
from SVclone.SVprocess import svp_load_data as load_data
from SVclone.SVprocess import count
from SVclone import load_data as svc_load
from SVclone import run_filter

blist = ''
bam = 'example_data/tumour_p80_DEL_sv_extract_sorted.bam'
snv_vcf = 'example_data/tumour_p80_DEL_snvs.vcf'
outdir = 'tumour_p80_DEL'
sample = 'tumour_p80_DEL'
svin_out = '%s/%s_svin.txt' % (outdir, sample)
svinfo_out = '%s/%s_svinfo.txt' % (outdir, sample)
cfg = 'svclone_config.ini'
pp_file = 'example_data/purity_ploidy.txt'
param_file = '%s/read_params.txt' % outdir
sv_filt_file = '%s/%s_filtered_svs.tsv' % (outdir, sample)
snv_filt_file = '%s/%s_filtered_snvs.tsv' % (outdir, sample)
svf_subset = '%s/%s_filtered_svs_subset.tsv' % (outdir, sample)
snvf_subset = '%s/%s_filtered_snvs_subset.tsv' % (outdir, sample)

if not os.path.exists(outdir):
    os.makedirs(outdir)

pi, ploidy = load_data.get_purity_ploidy(pp_file, sample, outdir)
rlen, insert, insert_std = load_data.get_read_params('', sample, outdir)
svs = load_data.load_input_simple('example_data/tumour_p80_DEL_svs_simple.txt', use_dir=False, class_field='')
sv_types = ['INV', 'DEL', 'DUP', 'INTDUP', 'TRX', 'INTRX']
threshold = 6

# cluster outputs
ass_prob_tbl = '%s_assignment_probability_table.txt' % sample
ccert        = '%s_cluster_certainty.txt' % sample
cn_out       = '%s_copynumber.txt' % sample
fit          = '%s_fit.txt' % sample
ml_cn        = '%s_most_likely_copynumbers.txt' % sample
mult         = '%s_multiplicity.txt' % sample
sc_str       = '%s_subclonal_structure.txt' % sample
n_clus       = '%s_number_of_clusters.txt' % sample
ccfs         = '%s_vaf_ccf.txt' % sample

Config = configparser.ConfigParser()
Config.read(cfg)

max_cn         = float(Config.get('BamParameters', 'max_cn'))
mean_cov       = float(Config.get('BamParameters', 'mean_cov'))
sc_len         = int(Config.get('SVcountParameters', 'sc_len'))
threshold      = int(Config.get('SVcountParameters', 'threshold'))
dna_gain_class = Config.get('SVclasses', 'dna_gain_class').split(',')
dna_loss_class = Config.get('SVclasses', 'dna_loss_class').split(',')
cdefs   = {'dna_gain_class': dna_gain_class, 'dna_loss_class': dna_loss_class}
min_dep = 8

consens_dtype = [('bp1_ca_right', int), ('bp1_ca_left', int), \
                 ('bp2_ca_right', int), ('bp2_ca_left', int)]

rlen = bamtools.estimateTagSize(bam)
inserts = bamtools.estimateInsertSizeDistribution(bam)
inserts = (max(rlen*2,inserts[0]),inserts[1])

max_ins = inserts[0]+(3*inserts[1]) #max fragment size = mean fragment len + (fragment std * 3)
max_dep = ((mean_cov*(max_ins*2))/rlen)*max_cn

clus_th   = {'percent': 0.01, 'absolute': 10}

class test(unittest.TestCase):

    def test_01_annotate_count(self):
        self.assertTrue(len(svs)==100)
        ca = np.zeros(len(svs), dtype=consens_dtype)
        sv_wdir, ca = annotate.infer_sv_dirs(svs, ca, bam, max_dep, sc_len, threshold, blist)
        self.assertTrue(np.all(sv_wdir['dir1'] == '+'))
        self.assertTrue(np.all(sv_wdir['dir2'] == '-'))

        sv_wdir = annotate.classify_svs(svs, threshold)
        self.assertTrue(np.all([svc in sv_types for svc in sv_wdir['classification']]))

        # write output for count steps
        annotate.write_svs(sv_wdir, svin_out)

    def test_02_count(self):
        rparams = count.get_params(cfg, bam, sample, outdir)
        split_reads, span_reads, anom_reads = count.extract_sv_info(svin_out, bam, rparams, svinfo_out)
        sv_df = svc_load.load_svs(svinfo_out)

        # test some simple properties to make sure data looks sensible
        self.assertTrue(np.all([s > 0 for s in sv_df.spanning.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.split1.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.split2.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.support.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.norm1.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.norm2.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.span_norm1.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.span_norm2.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.split_norm1.values]))
        self.assertTrue(np.all([s > 0 for s in sv_df.split_norm2.values]))
        self.assertTrue(np.all([s > 0 and s < 1 for s in sv_df.vaf1.values]))
        self.assertTrue(np.all([s > 0 and s < 1 for s in sv_df.vaf2.values]))

    def test_03_filter(self):
        self.assertTrue(pi > 0)
        self.assertTrue(pi <= 1)

        # TODO: add SVs to be filtered to example data, add dummy germline data
        sv_df = svc_load.load_svs(svinfo_out)
        sv_df_filt = run_filter.run_simple_filter(sv_df, rlen, inserts[0], 1, 1, -1, 8, False, [], '')
        self.assertTrue(len(sv_df_filt) == len(sv_df))
        sv_df = sv_df_filt

        # TODO: add copy-number testing/duplications adjustment testing
        maj_allele = round(ploidy/2) if round(ploidy) > 1 else 1
        min_allele = round(ploidy/2) if round(ploidy) > 1 else 0
        default_gtype = '%d,%d,1.0' % (maj_allele, min_allele)
        sv_df['gtype1'] = default_gtype
        sv_df['gtype2'] = default_gtype

        sv_df = run_filter.adjust_sv_read_counts(sv_df, pi, ploidy, min_dep, rlen, Config)
        sv_df.to_csv(sv_filt_file, sep='\t', index=False, na_rep='')

        sv_df_filt = pd.read_csv(sv_filt_file, delimiter='\t', dtype=None, header=0, low_memory=False)
        sv_df_filt = pd.DataFrame(sv_df).fillna('')
        self.assertTrue(len(sv_df) == len(sv_df_filt))

        snv_df = svc_load.load_snvs_mutect(snv_vcf, 'sample')
        snv_df_filt = run_filter.run_simple_snv_filter(snv_df, min_dep, blist, False, [])
        self.assertTrue(len(snv_df_filt) == len(snv_df))

        snv_df = snv_df_filt
        snv_df['gtype'] = default_gtype

        snv_df.to_csv(snv_filt_file, sep='\t', index=False, na_rep='')

if __name__ == '__main__':
    nose2.main()
