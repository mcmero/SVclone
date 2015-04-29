#!/usr/bin/env python

'''
Commandline input for running SV 
'''

import argparse
import phylo_sv
import ipdb

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument("-i","--input",dest="proc_svs",default="")
parser.add_argument("-b","--bam",dest="bam",default="")
parser.add_argument("-o","--out",dest="out")
args = parser.parse_args()

svs = args.svin
bam = args.bam
out = args.out

if __name__ == '__main__':
    phylo_sv.build_tree(svs,bam,out)

