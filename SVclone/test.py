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

blist = ''
bam = 'example_data/tumour_p80_DEL_sv_extract_sorted.bam'
outdir = 'tumour_p80_DEL'
sample = 'tumour_p80_DEL'
svin_out = '%s/%s_svin.txt' % (outdir, sample)
svinfo_out = '%s/%s_svinfo.txt' % (outdir, sample)
cfg = 'svclone_config.ini'

svs = load_data.load_input_simple('example_data/tumour_p80_DEL_svs_simple.txt', use_dir=False, class_field='')
sv_types = ['INV', 'DEL', 'DUP', 'INTDUP', 'TRX', 'INTRX']

Config = ConfigParser.ConfigParser()
Config.read(cfg)

max_cn   = int(Config.get('BamParameters', 'max_cn'))
mean_cov = int(Config.get('BamParameters', 'mean_cov'))
sc_len   = int(Config.get('SVcountParameters', 'sc_len'))
threshold = int(Config.get('SVcountParameters', 'threshold'))

consens_dtype = [('bp1_ca_right', int), ('bp1_ca_left', int), \
                 ('bp2_ca_right', int), ('bp2_ca_left', int)]

rlen = bamtools.estimateTagSize(bam)
inserts = bamtools.estimateInsertSizeDistribution(bam)
inserts = (max(rlen*2,inserts[0]),inserts[1])

max_ins = inserts[0]+(3*inserts[1]) #max fragment size = mean fragment len + (fragment std * 3)
max_dep = ((mean_cov*(max_ins*2))/rlen)*max_cn

if not os.path.exists(outdir):
    os.makedirs(outdir)

class test_identify(unittest.TestCase):

    def test_annotate_count(self):
        self.assertTrue(len(svs)==100)
        ca = np.zeros(len(svs), dtype=consens_dtype)
        sv_wdir, ca = identify.infer_sv_dirs(svs, ca, bam, max_dep, sc_len, threshold, blist)
        self.assertTrue(np.all(sv_wdir['dir1'] == '+'))
        self.assertTrue(np.all(sv_wdir['dir2'] == '-'))

        sv_wdir = identify.classify_svs(svs)
        self.assertTrue(np.all([svc in sv_types for svc in sv_wdir['classification']]))

        # write output for count steps
        identify.write_svs(sv_wdir, svin_out)

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

if __name__ == '__main__':
    nose2.main()
