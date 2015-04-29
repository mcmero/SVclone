#!/usr/bin/env python

'''
Commandline input for running SV post-processing script
'''

import argparse
import phylo_sv
import ipdb

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument("-i","--input",dest="svin",default="")
parser.add_argument("-b","--bam",dest="bam",default="")
parser.add_argument("-hd","--header",dest="header",default="")
parser.add_argument("-o","--out",dest="out")
parser.add_argument("-d","--depth",dest="mean_depth")
args = parser.parse_args()

svin = args.svin
bam = args.bam
out = args.out
hd_cfg = args.header
mean_dp = float(args.mean_depth)

if __name__ == '__main__':
    phylo_sv.proc_svs(svin,bam,out,hd_cfg,mean_dp)

