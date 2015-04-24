#!/usr/bin/env python

'''
Author: Marek Cmero
Commandline input for running sv_proc
'''

import argparse
import proc_svs
import ipdb

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument("svin")
parser.add_argument("bam")
parser.add_argument("header")
parser.add_argument("out")
parser.add_argument("mean_depth")
args = parser.parse_args()

svin = args.svin
bam = args.bam
out = args.out
hd_cfg = args.header
mean_dp = float(args.mean_depth)

if __name__ == '__main__':
    proc_svs.run(svin,bam,out,hd_cfg,mean_dp)
