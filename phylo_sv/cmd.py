#!/usr/bin/env python

'''
Commandline input for running SV 
'''

import argparse
import phylo_sv
import sys

parser = argparse.ArgumentParser(prefix_chars='--')
parser.add_argument("-s","--samples",dest="samples",
    help="Required: Sample names (comma separated if multiple), not including germline.")
parser.add_argument("-i","--input",dest="procd_svs",
    help="Required: Processed structural variation input (comma separated if multiple).")
parser.add_argument("-g","--germline",dest="germline",default="",
    help="Germline SVs. If not provided, will assume all SVs are somatic.")
parser.add_argument("-c","--cnvs",dest="cnvs",default="",
    help="Phased copy-number states from Battenberg (comma separated if multiple). " + \
            "If not provided, all SVs assumed copy-neutral.")
parser.add_argument("-r","--readlen",dest="rlen",default="100",
    help="Read length of samples. May be comma separated per sample. "+ \
            "If one value is given, read len is the same for all samples")
parser.add_argument("-v","--insert",dest="insert",default="",
    help="Mean insert size, where value = fragment size - (2 * read len). " +
            "Comma separated if multuple samples.")
parser.add_argument("-p","--purity",dest="pi",default="1.",
    help="Tumour purities for all samples given. A single parameter assumes " +
            "uniform purity for all samples. No parameter assumes 100% purity.")
parser.add_argument("-o","--out",dest="out",
    help="Base name for outputting phylogeny info.")
args = parser.parse_args()

samples = args.samples
svs     = args.procd_svs
gml     = args.germline
cnvs    = args.cnvs
out     = args.out
rlen    = args.rlen
insert  = args.insert
pi      = args.pi

def proc_arg(arg,n_args=1,of_type=str):
    arg = str.split(arg,',')
    arg = arg * n_args if len(arg)==1 else arg
    return map(of_type,arg)

if __name__ == '__main__':
    try:
        if insert=="":
            print("Inserts not provided, assuming insert length equals read length")
            args.insert=rlen
        samples = proc_arg(samples)
        n = len(samples)
        svs = proc_arg(svs)
        cnvs = proc_arg(cnvs)
        if len(svs)!=n or len(cnvs)!=n:
            raise ValueError
        rlen = proc_arg(rlen,n,int)
        insert = proc_arg(insert,n,float)
        pi = proc_arg(pi,n,float)
    except ValueError:
        print "Invalid arguments. Check arguments with -h or --help and try again."
        sys.exit
    phylo_sv.run(samples,svs,gml,cnvs,rlen,insert,pi,out)    
