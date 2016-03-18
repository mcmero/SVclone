#!/usr/bin/env python

'''
Commandline input for running SV
'''

from SVClone import run_filter
from SVClone import run_clus
from SVClone.SVProcess import identify
from SVClone.SVProcess import count

import argparse
import numpy as np

parser = argparse.ArgumentParser(prog='SVClone')

parser.add_argument('--version', action='version', version='SVClone-1.0.0')

subparsers = parser.add_subparsers()

##########################################################################################################

identify_parser = subparsers.add_parser('identify', help='Extract directions and SV classifications')

identify_parser.add_argument("-cfg","--config",dest="cfg",default="svclone_config.ini",
                    help="Config file.")

identify_parser.add_argument("-i","--input",dest="svin",required=True,
                    help="Structural variants input file. See README for input format")

identify_parser.add_argument("-b","--bam",dest="bam",required=True,
                    help="Corresponding indexed BAM file")

identify_parser.add_argument("-o","--out",dest="out",required=True,
                    help='''Output base name. May contain directories. Will create pre-processed output as 
                    <name>_svin.txt''')

#identify_parser.add_argument("-cn","--max_cn",dest="max_cn",default=10,type=int,
#                    help='''Optional: maximum expected copy-number. Will skip the processing of any areas 
#                    where reads > average coverage * max_cn''')
#
#identify_parser.add_argument("-d","--mean_depth",dest="mean_dp",type=float,default=50,
#                    help='''Average coverage for BAM file in covered region. May be calculated across 
#                    binned intervals and may be approximate''')

identify_parser.add_argument("--simple",dest="simple_svs",action="store_true",
                    help="Whether sv input is in a simple format (see README), otherwise VCF format is assumed.")

identify_parser.add_argument("--socrates",dest="socrates",action="store_true",
                    help="Whether sv input is 'Socrates' SV caller input.")

identify_parser.add_argument("--sv_class_field",dest="class_field",default="",
                    help="Use existing classification field, specify the field name")

identify_parser.add_argument("--use_dir",dest="use_dir",action="store_true",
                    help="Whether to use breakpoint direction in the input file (where it must be supplied).")

identify_parser.add_argument("--filter_repeats",dest="filt_repeats",default="",
                    help='''Comma-separated repeat types to filter, if found at both sides of the breakpoint). 
                    SOCRATES INPUT ONLY.''')

identify_parser.add_argument("--min_mapq",dest="min_mapq",default=0,type=float,
                    help='''Filter out SVs with lower average MAPQ than this value. SOCRATES INPUT ONLY (default 0)''')

identify_parser.add_argument("--trust_sc_pos",dest="trust_sc_pos",action="store_true",
                    help='''Use specified breaks without checking for differing soft-clip consensus position. 
                    Cannot be skipped if directionality must be inferred. If your SV caller offsets 
                    breaks due to micro-homology, e.g. Socrates, using this option is not recommended.''')

identify_parser.add_argument("-r","--read_len",dest="rlen",default=-1,type=int,
                    help="Read length. If not specified, will be inferred")

identify_parser.set_defaults(func=identify.preproc_svs)

##########################################################################################################

count_parser = subparsers.add_parser('count', help='Count reads from called structural variations')

count_parser.add_argument("-cfg","--config",dest="cfg",default="svclone_config.ini",
                    help="Config file.")

count_parser.add_argument("-i","--input",dest="svin",required=True,
                   help="Structural variants input file. See README for input format")

count_parser.add_argument("-b","--bam",dest="bam",required=True,
                    help="Corresponding indexed BAM file")

count_parser.add_argument("-o","--out",dest="out",required=True,
                    help='''Output base name. May contain directories. Will create processed output as 
                    <name>_svinfo.txt.''')

#count_parser.add_argument("-d","--mean_depth",dest="mean_dp",type=float,default=50,
#                    help='''Average coverage for BAM file in covered region. May be calculated across 
#                    binned intervals and may be approximate''')
#
#count_parser.add_argument("-sc","--softclip",dest="sc_len",default=10,type=int,
#                    help='''Optional: minimum number of basepairs by which reads spanning the break are 
#                    considered support the breakpoint. Also affects number of base-pairs a normal read 
#                    must overlap the break to be counted. Default = 10''')
#
#count_parser.add_argument("-cn","--max_cn",dest="max_cn",default=10,type=int,
#                    help='''Optional: maximum expected copy-number. Will skip the processing of any areas 
#                    where reads > average coverage * max_cn''')

count_parser.add_argument("-r","--read_len",dest="rlen",default=-1,type=int,
                    help="Read length. If not specified, will be inferred")

count_parser.add_argument("-v","--insert_mean",dest="insert_mean",default=-1.,type=float,
                    help='''Mean insert length (this is the fragment length or template length). 
                    If not specified, will be inferred''')

count_parser.add_argument("--insert_std",dest="insert_std",default=-1.,type=float,
                    help="Standard deviation of insert length. If not specified, will be inferred")

count_parser.add_argument("--write_anomalous",dest="write_anom",action="store_true",
                    help='''Whether to output (as a bam) and accurately recount anomalous reads. 
                    Useful for diagnosing issues with read counting.''')

count_parser.set_defaults(func=count.proc_svs)

##########################################################################################################

filter_parser = subparsers.add_parser('filter', help='Filter output from process step')

filter_parser.add_argument("-cfg","--config",dest="cfg",default="svclone_config.ini",
                    help="Config file. Default: svclone_config.ini")

filter_parser.add_argument("-s","--sample",dest="sample",required=True,
                    help='''Sample name.
                    WARNING: if clustering using mutect SNVs, the sample name must match the sample name 
                    in the vcf file.''')

filter_parser.add_argument("-i","--input",default="",dest="procd_svs",
                    help="Required: Processed structural variation input (comma separated if multiple).")

filter_parser.add_argument("-g","--germline",dest="germline",default="",
                    help='''Germline SVs in output format from process step. If not provided, will 
                    assume all SVs are somatic.''')

filter_parser.add_argument("-c","--cnvs",dest="cnvs",default="",
                    help='''Phased copy-number states from Battenberg. If not provided, all SVs assumed copy-neutral.''')

filter_parser.add_argument("--min_depth",dest="min_dep",type=float,default=4,
                    help='''Filter out any variants with total depth below this value (default = 4). Applies to
                    SVs and SNVs.''')

filter_parser.add_argument("--params",dest="params_file",default="",
                    help='''Parameters file from processing step. If not supplied, the default search path 
                    is <outdir>/<sample>_params.txt. If the file does not exist, a read length and mean
                    insert length of 100 will be selected.''')

filter_parser.add_argument("--neutral",dest="neutral",action="store_true",
                    help="Keep only copy-number neutral SVs.")

filter_parser.add_argument("--snvs",dest="snvs",default="",type=str,
                    help="SNVs in VCF format to (optionally) compare the clustering with SVs.")

filter_parser.add_argument("--snv_format",dest="snv_format",
                    choices=['sanger','mutect','mutect_callstats'],default="sanger",
                    help='''Supplied SNV VCF is in the following input format: sanger (default) or mutect,
                    mutect_callstats is also an option for the non-VCF mutect output.''')

filter_parser.add_argument("-o","--outdir",dest="outdir",default=".",
                    help="Output directory. Default: current directory")

filter_parser.add_argument("-p","--purity",dest="pi",default="1.",
                    help='''Tumour puritie for given sample. Default is 100%% purity.''')

filter_parser.add_argument("-y","--ploidy",dest="ploidy",default="2.0",
                    help="Tumour ploidy; default = 2 (diploid).")

filter_parser.add_argument("--minsplit",dest="minsplit",default=1,
                    help="Require at least N split reads to keep SV (default = 1).")

filter_parser.add_argument("--minspan",dest="minspan",default=1,
                    help="Require at least N spanning reads to keep SV (default = 1).")

filter_parser.add_argument("--sizefilter",dest="sizefilter",default=-1,type=int,
                    help='''Filter out SVs below this size. By default, SVs below read length * 2 + 
                    mean insert size are filtered out''')

filter_parser.add_argument("--filter_outliers",dest="filter_outliers",action="store_true",
                    help='''Filter out SVs with depth values that are considers outliers, based on the 
                    copy-number adjusted distribution of depths. Use with caution, as this may filter 
                    out real subclonal SV events.''')

filter_parser.add_argument("--valid_chroms",dest="valid_chrs",action="store_true",
                    help='''Filters out SVs on non-valid chroms (i.e. mapping to contigs on non-
                    standard chromosomes. Can be specified in the parameters.py file.''')

filter_parser.set_defaults(func=run_filter.run)

##########################################################################################################

cluster_parser = subparsers.add_parser('cluster', help='Run clustering step')

cluster_parser.add_argument("-cfg","--config",dest="cfg",default="svclone_config.ini",
                    help="Config file.")

cluster_parser.add_argument("-s","--sample",dest="sample",required=True,
                    help='''Sample name.''')

cluster_parser.add_argument("-o","--outdir",dest="outdir",default=".",
                    help="Output directory. Default: current directory")

cluster_parser.add_argument("-n","--n_runs",dest="n_runs",default=1,type=int,
                    help="Number of times to run whole rounds of sampling.")

cluster_parser.add_argument("-t","--n_iter",dest="n_iter",default=10000,type=int,
                    help="Number of MCMC iterations.")

cluster_parser.add_argument("--params",dest="params_file",default="",
                    help='''Parameters file from processing step. If not supplied, the default search path 
                    is <outdir>/<sample>_params.txt''')

cluster_parser.add_argument("--burn",dest="burn",default=0,type=int,
                    help="Burn-in for MCMC (default 0.)")

cluster_parser.add_argument("--thin",dest="thin",default=1,type=int,
                    help="Thinning parameter for MCMC (default 1.)")

cluster_parser.add_argument("--no_plot",dest="plot",action="store_false",default=True,
                    help="Do not create output plot")

#cluster_parser.add_argument("--beta",dest="beta",default="0.9,1/0.9,2",type=str,
#                    help='''Comma separated; first two values determine the shape and scale (1/rate) 
#                    parameters used in the Dirichlet Processes' gamma function. The third value is the 
#                    initial value. Default values: "0.9,1/0.9,2" (shape = 0.9, scale = 1/0.9, init = 2)''')

cluster_parser.add_argument("--merge",dest="merge_clusts",action="store_true",
                    help="Set to merge clusters.")

cluster_parser.add_argument("--map",dest="use_map",action="store_true",
                    help='''Use maximum a-posteriori (MAP) fitting (may significantly increase runtime).
                    Using MAP will give access to the run-fit metrics: the Akaike information criterion (AIC) and
                    the Bayesian Information Criterion (BIC).''')

cluster_parser.add_argument("--cocluster",dest="cocluster",action="store_true",
                    help="Whether to cluster SNVs and SVs together.")

cluster_parser.add_argument("--no_adjust",action="store_true",
                    help='''Do not use adjusted normal reads for duplications, or adjusted supporting reads 
                    for inversions''') 

cluster_parser.set_defaults(func=run_clus.run_clustering)

##########################################################################################################

args = parser.parse_args()
args.func(args)
