import nose2
import unittest
import numpy as np
import ConfigParser
import os
import pandas as pd
import subprocess
from unittest import TestCase
from SVprocess import bamtools
from SVprocess import annotate
from SVprocess import load_data
from SVprocess import count
from SVclone import load_data as svc_load
from SVclone import run_filter
from SVclone import run_clus
from SVclone import post_assign

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
n_clus       = 'number_of_clusters.txt'
ccfs         = '%s_vaf_ccf.txt' % sample

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

        # subset files for post-assign tests
        sv_df[:50].to_csv(svf_subset, sep='\t', index=False, na_rep='')
        snv_df[:250].to_csv(snvf_subset, sep='\t', index=False, na_rep='')

    def test_04_cluster(self):
        sv_df = pd.read_csv(sv_filt_file, delimiter='\t', dtype=None, header=0, low_memory=False)
        sv_df = pd.DataFrame(sv_df).fillna('')
        sample_params, cluster_params, output_params = \
            svc_load.get_params_cluster_step(sample, cfg, outdir, pp_file, param_file, False, False)

        snv_df = pd.read_csv(snv_filt_file, delimiter='\t', dtype=None, header=0, low_memory=False)
        snv_df = pd.DataFrame(snv_df).fillna('')

        cluster_params['n_iter'] = 110
        cluster_params['burn'] = 10
        cluster_params['use_map'] = False
        cluster_params['cocluster'] = False

        # separate clustering test
        subprocess.call(['rm', '-rf', '%s/run0' % outdir])
        run_clus.cluster_and_process(sv_df, snv_df, 0, outdir, sample_params,
                                     cluster_params, output_params, [1234])

        sv1 = pd.read_csv('%s/run0/%s' % (outdir, ass_prob_tbl),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv2 = pd.read_csv('%s/run0/%s' % (outdir, ccert), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv3 = pd.read_csv('%s/run0/%s' % (outdir, cn_out), delimiter='\t', dtype=None, header=0, low_memory=False)
        #sv4 = pd.read_csv('%s/run0/%s' % (outdir, fit), delimiter='\t', dtype=None, header=None, low_memory=False)
        sv5 = pd.read_csv('%s/run0/%s' % (outdir, ml_cn), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv6 = pd.read_csv('%s/run0/%s' % (outdir, mult), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv7 = pd.read_csv('%s/run0/%s' % (outdir, sc_str), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv8 = pd.read_csv('%s/run0/%s' % (outdir, n_clus), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv9 = pd.read_csv('%s/run0/%s' % (outdir, ccfs),
                          delimiter='\t', dtype=None, header=0, low_memory=False)

        snv9 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ccfs), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv1 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ass_prob_tbl),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv2 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ccert), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv3 = pd.read_csv('%s/run0/snvs/%s' % (outdir, cn_out), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv5 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ml_cn), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv6 = pd.read_csv('%s/run0/snvs/%s' % (outdir, mult), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv7 = pd.read_csv('%s/run0/snvs/%s' % (outdir, sc_str), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv8 = pd.read_csv('%s/run0/snvs/%s' % (outdir, n_clus), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv9 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ccfs), delimiter='\t', dtype=None, header=0, low_memory=False)

        self.assertTrue(len(sv1) == len(sv_df))
        self.assertTrue(len(sv2) == len(sv_df))
        self.assertTrue(len(sv3) == len(sv_df))
        #self.assertTrue(len(sv4) == 9)
        self.assertTrue(len(sv5) == len(sv_df))
        self.assertTrue(len(sv6) == len(sv_df))
        self.assertTrue(len(sv7.columns) == 4)
        self.assertTrue(len(sv8) == 1)
        self.assertTrue(len(sv9) == len(sv_df))
        self.assertTrue(len(snv1) == len(snv_df))
        self.assertTrue(len(snv2) == len(snv_df))
        self.assertTrue(len(snv3) == len(snv_df))
        self.assertTrue(len(snv5) == len(snv_df))
        self.assertTrue(len(snv6) == len(snv_df))
        self.assertTrue(len(snv7.columns) == 4)
        self.assertTrue(len(snv8) == 1)
        self.assertTrue(len(snv9) == len(snv_df))

        # coclustering test
        subprocess.call(['rm', '-rf', '%s/run0' % outdir])
        cluster_params['cocluster'] = True

        run_clus.cluster_and_process(sv_df, snv_df, 0, outdir, sample_params,
                                        cluster_params, output_params, [1234])

        snv1 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ass_prob_tbl),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv2 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ccert), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv3 = pd.read_csv('%s/run0/snvs/%s' % (outdir, cn_out), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv5 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ml_cn), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv6 = pd.read_csv('%s/run0/snvs/%s' % (outdir, mult), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv7 = pd.read_csv('%s/run0/snvs/%s' % (outdir, sc_str), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv8 = pd.read_csv('%s/run0/snvs/%s' % (outdir, n_clus), delimiter='\t', dtype=None, header=0, low_memory=False)
        snv9 = pd.read_csv('%s/run0/snvs/%s' % (outdir, ccfs), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv1 = pd.read_csv('%s/run0/%s' % (outdir, ass_prob_tbl), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv2 = pd.read_csv('%s/run0/%s' % (outdir, ccert), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv3 = pd.read_csv('%s/run0/%s' % (outdir, cn_out), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv5 = pd.read_csv('%s/run0/%s' % (outdir, ml_cn), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv6 = pd.read_csv('%s/run0/%s' % (outdir, mult), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv7 = pd.read_csv('%s/run0/%s' % (outdir, sc_str), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv8 = pd.read_csv('%s/run0/%s' % (outdir, n_clus), delimiter='\t', dtype=None, header=0, low_memory=False)
        sv9 = pd.read_csv('%s/run0/%s' % (outdir, ccfs), delimiter='\t', dtype=None, header=0, low_memory=False)

        self.assertTrue(len(sv1) == len(sv_df))
        self.assertTrue(len(sv2) == len(sv_df))
        self.assertTrue(len(sv3) == len(sv_df))
        self.assertTrue(len(sv5) == len(sv_df))
        self.assertTrue(len(sv6) == len(sv_df))
        self.assertTrue(len(sv8) == 1)
        self.assertTrue(len(sv9) == len(sv_df))
        self.assertTrue(len(snv1) == len(snv_df))
        self.assertTrue(len(snv2) == len(snv_df))
        self.assertTrue(len(snv3) == len(snv_df))
        self.assertTrue(len(snv5) == len(snv_df))
        self.assertTrue(len(snv6) == len(snv_df))
        self.assertTrue(len(snv8) == 1)
        self.assertTrue(len(snv9) == len(snv_df))

        # TODO - add tests for run selection

    def test_05_post_assign(self):
        sv_df = pd.read_csv(sv_filt_file, delimiter='\t', dtype=None, header=0, low_memory=False)
        sv_df = pd.DataFrame(sv_df).fillna('')

        snv_df = pd.read_csv(snv_filt_file, delimiter='\t', dtype=None, header=0, low_memory=False)
        snv_df = pd.DataFrame(snv_df).fillna('')

        sv_filt_df = sv_df.head(50).copy()
        snv_filt_df = snv_df.head(250).copy()

        # rerun clustering with subset
        subprocess.call(['rm', '-rf', '%s/run0' % outdir])
        sample_params, cluster_params, output_params = \
            svc_load.get_params_cluster_step(sample, cfg, outdir, pp_file, param_file, False, False)

        cluster_params['n_iter'] = 100
        cluster_params['burn'] = 0
        cluster_params['cocluster'] = True
        cluster_params['use_map'] = False

        run_clus.cluster_and_process(sv_filt_df, snv_filt_df, 0, outdir, sample_params,
                                        cluster_params, output_params, [1234])

        rundir = '%s/run0' % outdir
        sv_to_assign = post_assign.get_var_to_assign(sv_df, sv_filt_df)

        self.assertTrue(len(sv_to_assign) == (len(sv_df) - len(sv_filt_df)))

        sv_to_assign = run_filter.adjust_sv_read_counts(sv_to_assign, pi, ploidy, 0, rlen, Config)
        post_assign.post_assign_vars(sv_to_assign, sv_filt_df, rundir, sample, sample_params, cluster_params)

        sv1 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, ass_prob_tbl),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv2 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, ccert),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv3 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, cn_out),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv5 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, ml_cn),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv6 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, mult),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv7 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, sc_str),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv8 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, n_clus),
                          delimiter='\t', dtype=None, header=0, low_memory=False)
        sv9 = pd.read_csv('%s/run0_post_assign/%s' % (outdir, ccfs),
                          delimiter='\t', dtype=None, header=0, low_memory=False)

        self.assertTrue(len(sv1) == len(sv_df))
        self.assertTrue(len(sv2) == len(sv_df))
        self.assertTrue(len(sv3) == len(sv_df))
        self.assertTrue(len(sv5) == len(sv_df))
        self.assertTrue(len(sv6) == len(sv_df))
        self.assertTrue(len(sv8) == 1)
        self.assertTrue(len(sv9) == len(sv_df))

        snv_to_assign = post_assign.get_var_to_assign(snv_df, snv_filt_df, snvs = True)
        post_assign.post_assign_vars(snv_to_assign, snv_filt_df, rundir, sample, sample_params, cluster_params, snvs = True)

        snv1 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, ass_prob_tbl),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv2 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, ccert),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv3 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, cn_out),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv5 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, ml_cn),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv6 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, mult),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv7 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, sc_str),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv8 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, n_clus),
                           delimiter='\t', dtype=None, header=0, low_memory=False)
        snv9 = pd.read_csv('%s/run0_post_assign/snvs/%s' % (outdir, ccfs),
                           delimiter='\t', dtype=None, header=0, low_memory=False)

        self.assertTrue(len(snv1) == len(snv_df))
        self.assertTrue(len(snv2) == len(snv_df))
        self.assertTrue(len(snv3) == len(snv_df))
        self.assertTrue(len(snv5) == len(snv_df))
        self.assertTrue(len(snv6) == len(snv_df))
        self.assertTrue(len(snv8) == 1)
        self.assertTrue(len(snv9) == len(snv_df))

    # TODO: add test for map/picking best run
    # TODO: add tests for cluster merging

if __name__ == '__main__':
    nose2.main()
