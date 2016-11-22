from __future__ import print_function

import os
import pandas as pd
import numpy as np
import ConfigParser

from . import load_data
from . import run_filter as filt
from SVprocess import load_data as svp_load

def get_sv_to_assign(sv_df, sv_filt_df):
    sv_to_assign = sv_df

    if len(sv_filt_df) > 0:
        ids = zip(sv_df.chr1.values, sv_df.pos1.values, sv_df.dir1,
                       sv_df.chr2.values, sv_df.pos2.values, sv_df.dir2)
        ids = [':'.join([str(y) for y in x]) for x in ids]

        filt_ids = zip(sv_filt_df.chr1.values, sv_filt_df.pos1.values, sv_filt_df.dir1,
                       sv_filt_df.chr2.values, sv_filt_df.pos2.values, sv_filt_df.dir2)
        filt_ids = [':'.join([str(y) for y in x]) for x in filt_ids]

        sv_to_assign = sv_df[[id not in filt_ids for id in ids]]

    return sv_to_assign

def run_post_assign(args):

    sample          = args.sample
    cfg             = args.cfg
    out             = sample if args.out == "" else args.out
    snv_file        = args.snv_file
    snv_format      = args.snv_format
    sv_file         = args.sv_file
    cnv_file        = args.cnv_file
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

    if out != '' and not os.path.exists(out):
        raise ValueError('Specified output directory does not exist!')

    if snv_file == '' and sv_file == '':
        raise ValueError('No variants specified!')

    if snv_file != '' and not os.path.exists(snv_file):
        raise OSError('Specified SNV file does not exist!')

    if sv_file != '' and not os.path.exists(sv_file):
        raise OSError('Specified SV file does not exist!')

    sv_filt_file = '%s/%s_filtered_svs.tsv' % (out, sample)
    snv_filt_file = '%s/%s_filtered_snvs.tsv' % (out, sample)

    sv_filt_df  = pd.DataFrame()
    snv_filt_df = pd.DataFrame()

    if os.path.exists(sv_filt_file):
        sv_filt_df = pd.read_csv(sv_filt_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        sv_filt_df = pd.DataFrame(sv_filt_df).fillna('')

    if os.path.exists(snv_filt_file):
        snv_filt_df = pd.read_csnv(snv_filt_file,delimiter='\t',dtype=None,header=0,low_memory=False)
        snv_filt_df = pd.DataFrame(snv_filt_df).fillna('')

    if len(sv_filt_df) == 0 and len(snv_filt_df) == 0:
        raise ValueError('Output directory filtered variant files do not exist or are empty!')

    param_file = '%s/read_params.txt' % out
    pp_file = '%s/purity_ploidy.txt' % out

    sample_params, cluster_params, output_params = \
                   load_data.get_params_cluster_step(sample, cfg, out, pp_file, param_file, XX, XY)
    rlen, insert, insert_std = svp_load.get_read_params(param_file, sample, out)
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
            sv_df = filter_germline(gml,sv_df,rlen,insert,gl_th)

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
        sv_to_assign = get_sv_to_assign(sv_df, sv_filt_df)

    if len(snv_df) > 0:
        snv_to_assign = get_snv_to_assign(snv_df, snv_filt_df)

