'''
Wrapper for post-assign R script
'''

import subprocess
import os

def run_postassign(args):
    sample          = args.sample
    out             = '%s/ccube_out/post_assign' % sample if args.out == "" else args.out
    joint           = '--joint' if args.joint else ''
    sv_rdata        = args.sv_rdata if args.sv_rdata != '' else '%s/ccube_out/%s_ccube_sv_results.RData' % (sample, sample)
    snv_rdata       = args.snv_rdata if args.snv_rdata != '' else '%s/ccube_out/snvs/%s_ccube_snv_results.RData' % (sample, sample)

    if out!='' and not os.path.exists(out):
        os.makedirs(out)

    if not os.path.exists(sv_rdata) or not os.path.exists(snv_rdata):
        raise OSError('One of SNV and/or SV RData input files not found!')
    
    dirname = os.path.dirname(os.path.abspath(__file__))
    subprocess.call(['Rscript', '%s/post_assign.R' % dirname,
                     sv_rdata, snv_rdata, out, sample, joint])
