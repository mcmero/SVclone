import nose2
import unittest
import numpy as np
import ConfigParser
import os
import pandas as pd
from unittest import TestCase
from SVprocess import bamtools
from SVprocess import identify
from SVprocess import load_data
from SVprocess import count
from SVclone import load_data as svc_load
from SVclone import run_filter
from SVclone import run_clus

blist = ''
bam = 'example_data/tumour_p80_DEL_sv_extract_sorted.bam'
outdir = 'tumour_p80_DEL'
sample = 'tumour_p80_DEL'
svin_out = '%s/%s_svin.txt' % (outdir, sample)
svinfo_out = '%s/%s_svinfo.txt' % (outdir, sample)
cfg = 'svclone_config.ini'
pp_file = 'example_data/purity_ploidy.txt'
param_file = '%s/read_params.txt' % outdir
sv_filt_file = '%s/%s_filtered_svs.tsv' % (outdir, sample)

if not os.path.exists(outdir):
    os.makedirs(outdir)

pi, ploidy = load_data.get_purity_ploidy(pp_file, sample, outdir)
rlen, insert, insert_std = load_data.get_read_params('', sample, outdir)
svs = load_data.load_input_simple('example_data/tumour_p80_DEL_svs_simple.txt', use_dir=False, class_field='')
sv_types = ['INV', 'DEL', 'DUP', 'INTDUP', 'TRX', 'INTRX']

# cluster outputs
ass_prob_tbl = '%s/run0/%s_assignment_probability_table.txt' % (outdir, sample)
ccert        = '%s/run0/%s_cluster_certainty.txt' % (outdir, sample)
cn_out       = '%s/run0/%s_copynumber.txt' % (outdir, sample)
fit          = '%s/run0/%s_fit.txt' % (outdir, sample)
ml_cn        = '%s/run0/%s_most_likely_copynumbers.txt' % (outdir, sample)
mult         = '%s/run0/%s_multiplicity.txt' % (outdir, sample)
sc_str       = '%s/run0/%s_subclonal_structure.txt' % (outdir, sample)
n_clus       = '%s/run0/number_of_clusters.txt' % outdir

Config = ConfigParser.ConfigParser()
Config.read(cfg)

max_cn   = int(Config.get('BamParameters', 'max_cn'))
mean_cov = int(Config.get('BamParameters', 'mean_cov'))
sc_len   = int(Config.get('SVcountParameters', 'sc_len'))
threshold = int(Config.get('SVcountParameters', 'threshold'))
min_dep = 8

consens_dtype = [('bp1_ca_right', int), ('bp1_ca_left', int), \
                 ('bp2_ca_right', int), ('bp2_ca_left', int)]

rlen = bamtools.estimateTagSize(bam)
inserts = bamtools.estimateInsertSizeDistribution(bam)
inserts = (max(rlen*2,inserts[0]),inserts[1])

max_ins = inserts[0]+(3*inserts[1]) #max fragment size = mean fragment len + (fragment std * 3)
max_dep = ((mean_cov*(max_ins*2))/rlen)*max_cn

class test_identify(unittest.TestCase):

    def test_01_annotate_count(self):
        self.assertTrue(len(svs)==100)
        ca = np.zeros(len(svs), dtype=consens_dtype)
        sv_wdir, ca = identify.infer_sv_dirs(svs, ca, bam, max_dep, sc_len, threshold, blist)
        self.assertTrue(np.all(sv_wdir['dir1'] == '+'))
        self.assertTrue(np.all(sv_wdir['dir2'] == '-'))

        sv_wdir = identify.classify_svs(svs)
        self.assertTrue(np.all([svc in sv_types for svc in sv_wdir['classification']]))

        # write output for count steps
        identify.write_svs(sv_wdir, svin_out)

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

        # TODO: add SNV loading test

    def test_04_cluster(self):
        sv_df = pd.read_csv(sv_filt_file, delimiter='\t', dtype=None, header=0, low_memory=False)
        sv_df = pd.DataFrame(sv_df).fillna('')
        sample_params, cluster_params, output_params = \
            svc_load.get_params_cluster_step(sample, cfg, outdir, pp_file, param_file, False, False)
        run_clus.cluster_and_process(sv_df, pd.DataFrame(), 0, outdir, sample_params,
                                        cluster_params, output_params, [1234])

        f1 = pd.read_csv(ass_prob_tbl, delimiter='\t', dtype=None, header=0, low_memory=False)
        f2 = pd.read_csv(ccert, delimiter='\t', dtype=None, header=0, low_memory=False)
        f3 = pd.read_csv(cn_out, delimiter='\t', dtype=None, header=0, low_memory=False)
        f4 = pd.read_csv(fit, delimiter='\t', dtype=None, header=None, low_memory=False)
        f5 = pd.read_csv(ml_cn, delimiter='\t', dtype=None, header=0, low_memory=False)
        f6 = pd.read_csv(mult, delimiter='\t', dtype=None, header=0, low_memory=False)
        f7 = pd.read_csv(sc_str, delimiter='\t', dtype=None, header=0, low_memory=False)
        f8 = pd.read_csv(n_clus, delimiter='\t', dtype=None, header=0, low_memory=False)

        self.assertTrue(len(f1) == len(sv_df))
        self.assertTrue(len(f2) == len(sv_df))
        self.assertTrue(len(f3) == len(sv_df))
        self.assertTrue(len(f4) == 9)
        self.assertTrue(len(f5) == len(sv_df))
        self.assertTrue(len(f6) == len(sv_df))
        self.assertTrue(len(f7.columns) == 4)
        self.assertTrue(len(f8) == 1)

        # TODO: more thorough tests checking output
        # TODO: add tests for cluster merging, coclustering, sv simulation

if __name__ == '__main__':
    nose2.main()
