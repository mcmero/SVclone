#!/usr/bin/env python

'''
Commandline input for running SV post-processing script
'''

import argparse
import numpy as np
from . import preprocess

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument("-i","--input",dest="svin",
        help="Structural variants input file. See README for input format")
parser.add_argument("-b","--bam",dest="bam",
        help="Corresponding indexed BAM file")
parser.add_argument("-hd","--header",dest="header",default="",
        help="Config file specifying header columns. If not provided, default column names assumed (see README)")
parser.add_argument("-o","--out",dest="out",
        help="Output base name. Will create processed output as <name>.txt and database output as <name>.db")
parser.add_argument("-d","--depth",dest="mean_dp",type=float,
        help="Average coverage for BAM file in covered region. May be calculated across binned intervals.")
parser.add_argument("-sc","--softclip_bp",dest="sc_len",default=25,type=int,
        help="Optional: minimum number of basepairs by which reads spanning the break are considered support " + \
             "the breakpoint. Also affects number of base-pairs a normal read must overlap the break to be " + \
             "counted. Default = 25")
parser.add_argument("-cn","--max_cn",dest="max_cn",default=15,type=int,
        help="Optional: maximum expected copy-number. Will skip the processing of any areas where" + \
             "# reads > average coverage * max_cn")
args = parser.parse_args()

svin    = args.svin
bam     = args.bam
out     = args.out
hd_cfg  = args.header
mean_dp = float(args.mean_dp)
sc_len  = int(args.sc_len)
max_cn  = int(args.max_cn)

if __name__ == '__main__':
    preprocess.proc_svs(svin,bam,out,hd_cfg,mean_dp,sc_len,max_cn)

