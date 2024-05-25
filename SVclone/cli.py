#!/usr/bin/env python

'''
Commandline input for running SVclone
'''

from SVclone import run_filter
from SVclone import run_clus
from SVclone import run_postassign
from SVclone.SVprocess import annotate
from SVclone.SVprocess import count

import argparse
import numpy as np
import pkg_resources

def main():

	parser = argparse.ArgumentParser(prog='SVclone')

	parser.add_argument('--version', action='version', version=pkg_resources.get_distribution("SVclone").version)

	subparsers = parser.add_subparsers()

	##########################################################################################################

	annotate_parser = subparsers.add_parser('annotate', help='Extract directions and SV classifications')

	annotate_parser.add_argument("-cfg","--config",dest="cfg",default="svclone_config.ini",
						help="Config file.")

	annotate_parser.add_argument("-i","--input",dest="svin",required=True,
						help="Structural variants input file. See README for input format")

	annotate_parser.add_argument("-b","--bam",dest="bam",required=True,
						help="Corresponding indexed BAM file")

	annotate_parser.add_argument("-s","--sample",dest="sample",required=True,
						help='''Sample name. Output is written to <out_dir>/<sample>_svin.txt.''')

	annotate_parser.add_argument("-o","--out",dest="out",default="",
						help='''Output directory. Sample name by default.''')

	annotate_parser.add_argument("--sv_format",dest="sv_format",choices=['vcf','simple','socrates'],default='vcf',
						help="Possible SV input formats: vcf, simple, socrates")

	annotate_parser.add_argument("--blacklist", dest="blist", default="",
						help='''Takes a file in BED format as an argument. Skip processing of any break-pairs
						where either SV break-end overlaps an interval specified in the supplied bed file.''')

	annotate_parser.set_defaults(func=annotate.preproc_svs)

	##########################################################################################################

	count_parser = subparsers.add_parser('count', help='Count reads from called structural variations')

	count_parser.add_argument("-cfg","--config",dest="cfg",default="svclone_config.ini",
						help="Config file.")

	count_parser.add_argument("-i","--input",dest="svin",required=True,
					   help="Structural variants input file. See README for input format")

	count_parser.add_argument("-b","--bam",dest="bam",required=True,
						help="Corresponding indexed BAM file")

	count_parser.add_argument("-s","--sample",dest="sample",required=True,
						help='''Sample name. Output is written to <out_dir>/<sample>_svinfo.txt.''')

	count_parser.add_argument("-o","--out",dest="out",default="",
						help='''Output directory. Default: sample name.''')

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

	filter_parser.add_argument("--params",dest="param_file",default="",
						help='''Read parameters file containing read length, average insert and insert standard
						deviation (see README). If not supplied, the default search path is <outdir>/<sample>_params.txt.
						If the file does not exist, a read length of 100 and mean insert length of 300 will be allocated.''')

	filter_parser.add_argument("--snvs",dest="snvs",default="",type=str,
						help="SNVs in VCF format to (optionally) compare the clustering with SVs.")

	filter_parser.add_argument("--snv_format",dest="snv_format",
						choices=['sanger','mutect','mutect_callstats','consensus','multisnv'],default="sanger",
						help='''Supplied SNV VCF is in the following input format: sanger (default), mutect,
						consensus (PCAWG) or mutect_callstats (non-VCF).''')

	filter_parser.add_argument("-o","--out",dest="out",default="",
						help='''Output directory. Default: sample name.''')

	filter_parser.add_argument("-p","--purity_ploidy",dest="pp_file",default="",
						help='''Tumour purity ploidy file. See README for format. The default file path is
						<outdir>/purity_ploidy.txt. If not found, default purity = 1 (100%%); default ploidy = 2.''')

	filter_parser.add_argument("--blacklist", dest="blist", default="",
						help='''Takes a file in BED format as an argument. Filter out any break-pairs where
						either SV break-end overlaps an interval specified in the supplied bed file.''')

	filter_parser.set_defaults(func=run_filter.run)

	##########################################################################################################

	cluster_parser = subparsers.add_parser('cluster', help='Run clustering step')

	cluster_parser.add_argument("-cfg","--config",dest="cfg",default="svclone_config.ini",
						help="Config file.")

	cluster_parser.add_argument("-i","--input",default="",dest="sv_file",
						help="Filtered structural variant input from filter step. Default loc: <outdir>/<sample>_filtered_svs.tsv")

	cluster_parser.add_argument("-s","--sample",dest="sample",required=True,
						help='''Sample name.''')

	cluster_parser.add_argument("-o","--out",dest="out",default="",
						help="Output directory. Default: sample name.")

	cluster_parser.add_argument("--params",dest="param_file",default="",
						help='''Read parameters file containing read length, average insert and insert standard
						deviation (see README). If not supplied, the default search path is <outdir>/<sample>_params.txt.
						If the file does not exist, a read length of 100 and mean insert length of 300 will be allocated.''')

	cluster_parser.add_argument("-p","--purity_ploidy",dest="pp_file",default="",
						help='''Tumour purity ploidy file. See README for format. The default file path is
						<outdir>/purity_ploidy.txt. If not found, default purity = 1 (100%%); default ploidy = 2.''')

	cluster_parser.add_argument("--snvs",dest="snv_file",default="",
						help="To specify filtered SNVs output from Filter Step. Default loc: <outdir>/<sample>_filtered_snvs.tsv")

	cluster_parser.add_argument("--subsample",dest="subsample",default=0,type=int,
						help="Subsample N SNVs from total filtered output to use for clustering.")

	cluster_parser.add_argument("--ss_seeds", dest="ss_seeds", default="",
						help='''Integer seeds to set seeds for replicability of subsampling.''')

	cluster_parser.add_argument("--XX",dest="XX",action="store_true",
						help="Specify XX genotype. (Overwrites config file, sets male to False.)")

	cluster_parser.add_argument("--XY",dest="XY",action="store_true",
						help="Specify XY genotype. (Overwrites config file, sets male to True.)")

	cluster_parser.set_defaults(func=run_clus.run_clustering)

	##########################################################################################################

	postassign_parser = subparsers.add_parser('postassign', help='Run post-assign step')

	postassign_parser.add_argument("-s","--sample",dest="sample",required=True,
						help='''Sample name.''')

	postassign_parser.add_argument("-o","--out",dest="out",default="",
						help="Output directory. Default: <sample>/ccube_out/post_assign.")

	postassign_parser.add_argument("-j","--joint",dest="joint",action="store_true",
						help="Indicates that variants will be post-assigned to a joint SV + SNV model")

	postassign_parser.add_argument("--svs",dest="sv_rdata",default="",
						help="SV RData file from postassign step output.")

	postassign_parser.add_argument("--snvs",dest="snv_rdata",default="",
						help="SNV RData file from postassign step output.")

	postassign_parser.set_defaults(func=run_postassign.run_postassign)

	##########################################################################################################

	args = parser.parse_args()
	args.func(args)




if __name__ == '__main__':
	main()
