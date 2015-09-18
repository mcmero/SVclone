#!/usr/bin/env python

'''
Commandline input for running SV
'''

from SVClone import run_filter
from SVClone import run_clus

import argparse
import numpy as np

parser = argparse.ArgumentParser(prog='SVClone')

parser.add_argument('--version', action='version', version='SVClone-0.1.0')

subparsers = parser.add_subparsers()

##########################################################################################################
filter_parser = subparsers.add_parser('filter', help='Filter output from process step')

filter_parser.add_argument("-s","--samples",dest="sample",
                    help='''Required: Sample name (comma separated if multiple), not including germline.
                    WARNING: if clustering using mutect SNVs, the sample name must match the sample name 
                    in the vcf file.''')

filter_parser.add_argument("-i","--input",default="",dest="procd_svs",
                    help="Required: Processed structural variation input (comma separated if multiple).")

filter_parser.add_argument("-g","--germline",dest="germline",default="",
                    help="Germline SVs. If not provided, will assume all SVs are somatic.")

filter_parser.add_argument("-c","--cnvs",dest="cnvs",default="",
                    help='''Phased copy-number states from Battenberg (comma separated if multiple).If not 
                    provided, all SVs assumed copy-neutral.''')

filter_parser.add_argument("-r","--readlen",dest="rlen",default="100",
                    help='''Read length of samples. May be comma separated per sample. If one value is given, 
                    read len is the same for all samples''')

filter_parser.add_argument("-v","--insert",dest="insert",default="",
                    help='''Mean insert size, where value = fragment size - (2 * read len).
                    Comma separated if multuple samples.''')

filter_parser.add_argument("--neutral",dest="neutral",action="store_true",
                    help="Keep only copy-number neutral SVs.")

filter_parser.add_argument("--snvs",dest="snvs",default="",type=str,
                    help="SNVs in VCF format to (optionally) compare the clustering with SVs.")

filter_parser.add_argument("--snv_format",dest="snv_format",
                    choices=['sanger','mutect','mutect_callstats'],default="sanger",
                    help='''Supplied SNV VCF is in the following input format: sanger (default), mutect or 
                    mutect_callstats.''')

filter_parser.add_argument("-o","--outdir",dest="outdir",default=".",
                    help="Output directory. Default: current directory")

filter_parser.add_argument("-p","--purity",dest="pi",default="1.",
                    help='''Tumour purities for all samples given. A single parameter assumes
                    uniform purity for all samples. No parameter assumes 100% purity.''')

filter_parser.add_argument("-y","--ploidy",dest="ploidy",default="2.0",
                    help="Tumour ploidy; default = 2 (diploid).")

filter_parser.set_defaults(func=run_filter.run)

##########################################################################################################

cluster_parser = subparsers.add_parser('cluster', help='Run clustering step')

cluster_parser.add_argument("-s","--samples",dest="sample",
                    help='''Required: Sample name (comma separated if multiple), not including germline.
                    WARNING: if clustering using mutect SNVs, the sample name must match the sample name 
                    in the vcf file.''')

cluster_parser.add_argument("-o","--outdir",dest="outdir",default=".",
                    help="Output directory. Default: current directory")

cluster_parser.add_argument("-n","--n_runs",dest="n_runs",default=1,type=int,
                    help="Number of times to run whole rounds of sampling.")

cluster_parser.add_argument("-t","--n_iter",dest="n_iter",default=10000,type=int,
                    help="Number of MCMC iterations.")

cluster_parser.add_argument("--burn",dest="burn",default=0,type=int,
                    help="Burn-in for MCMC (default 0.)")

cluster_parser.add_argument("--thin",dest="thin",default=1,type=int,
                    help="Thinning parameter for MCMC (default 1.)")

cluster_parser.add_argument("--plot",dest="plot",action="store_true",
                    help="Plot traces and clusters.")

cluster_parser.add_argument("--beta",dest="beta",default="8,1/0.05,0.1",type=str,
                    help='''Comma separated; first two values etermine the parameters used for
                    Dirichlet Processes' gamma function. Third value determines the starting value.''')

cluster_parser.add_argument("--merge",dest="merge_clusts",action="store_true",
                    help="Set to merge clusters.")

cluster_parser.add_argument("--map",dest="use_map",action="store_true",
                    help="Use maximum a-posteriori fitting (may significantly increase runtime).")

cluster_parser.add_argument("--cocluster",dest="cocluster",action="store_true",
                    help="Whether to cluster SNVs and SVs together.")

cluster_parser.set_defaults(func=run_clus.run_clustering)

##########################################################################################################

args = parser.parse_args()
args.func(args)
