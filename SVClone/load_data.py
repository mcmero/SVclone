import pandas as pd
import numpy as np
import vcf

def load_svs(sv_file):
    dat = pd.read_csv(sv_file,delimiter='\t',dtype=None, low_memory=False)
    sv_df = pd.DataFrame(dat)
    sv_df['norm_mean'] = map(np.mean,zip(sv_df['norm1'].values,sv_df['norm2'].values))
    return sv_df

def load_cnvs(cnv):
    cnv = pd.read_csv(cnv,delimiter='\t',dtype=None)
    cnv_df = pd.DataFrame(cnv)
    cnv_df['chr'] = map(str,cnv_df['chr'])
    
    try:
        gtypes = cnv_df['nMaj1_A'].map(str) + ',' + \
                 cnv_df['nMin1_A'].map(str) + ',' + \
                 cnv_df['frac1_A'].map(str)

        # join subclonal genotypes
        subclonal = cnv_df['frac1_A']!=1
        cnv_sc    = cnv_df[subclonal]
        gtypes[subclonal] = gtypes[subclonal] + '|' + \
                            cnv_sc['nMaj2_A'].map(str) + ',' + \
                            cnv_sc['nMin2_A'].map(str) + ',' + \
                            cnv_sc['frac2_A'].map(str)
        
        cnv_df['gtype'] = gtypes
        select_cols = ['chr','startpos','endpos','gtype']
        return cnv_df[select_cols]
    except KeyError:
        raise Exception('CNV file column names not recognised. Is the input a Battenberg output file?')

def load_snvs_mutect_callstats(snvs):
    snv_df = pd.DataFrame(pd.read_csv(snvs,delimiter='\t',low_memory=False,dtype=None,comment="#"))
    snv_df = snv_df[snv_df['judgement']=='KEEP']   
    
    snv_out = {'chrom' : snv_df.contig.map(str),
               'pos'   : snv_df.position.map(int),
               'gtype' : '',
               'ref'   : snv_df.t_ref_sum.map(float),
               'var'   : snv_df.t_alt_sum.map(float)}

    snv_out = pd.DataFrame(snv_out)
    snv_out = snv_out[['chrom','pos','gtype','ref','var']]
    return snv_out
    
def load_snvs_mutect(snvs,sample):
    vcf_reader = vcf.Reader(filename=snvs)
    snv_dtype = [('chrom','S50'),('pos',int),('gtype','S50'),('ref',float),('var',float)]
    snv_df = np.empty([0,5],dtype=snv_dtype)

    samples = vcf_reader.samples
    if len(samples)==0:
        raise Exception('No samples found in VCF!')
    elif not np.any(np.array(samples)==sample):
        print('Warning, sample not found in VCF, selecting first sample')
        sample = samples[0]

    for record in vcf_reader:
        if record.FILTER is not None:
            if len(record.FILTER)>0:
                continue
        ad = record.genotype(sample)['AD']
        ref_reads, variant_reads = float(ad[0]), float(ad[1])
        total_reads = ref_reads + variant_reads
        if variant_reads!=0:
            tmp = np.array((record.CHROM,record.POS,'',ref_reads,variant_reads),dtype=snv_dtype)
            snv_df = np.append(snv_df,tmp)

    return pd.DataFrame(snv_df)

def load_snvs_sanger(snvs):
    vcf_reader = vcf.Reader(filename=snvs)
    snv_dtype = [('chrom','S50'),('pos',int),('gtype','S50'),('ref',float),('var',float)]
    snv_df = np.empty([0,5],dtype=snv_dtype)

    #code adapted from: https://github.com/morrislab/phylowgs/blob/master/parser/create_phylowgs_inputs.py
    samples = vcf_reader.samples
    if samples[0]!='NORMAL' or samples[1]!='TUMOUR':
        raise Exception('VCF SNV file is of invalid format. Expected "NORMAL" and "TUMOUR" samples.')

    for record in vcf_reader:
        # get most likely genotypes
        genotypes = [record.INFO['TG'], record.INFO['SG']]        
        if len(genotypes)==0:
            continue
        if record.FILTER is not None:
            if len(record.FILTER)>0:
                continue
    
        variant_set = set()
        reference_nt = ''
        while len(variant_set) == 0:
            if len(genotypes) == 0:
                break
                #raise Exception('No more genotypes to find variant_nt in for %s' % variant)
            gt = genotypes.pop(0)
            normal_gt, tumour_gt = gt.split('/')
            if normal_gt[0] == normal_gt[1]:
                reference_nt = normal_gt[0]
                variant_set = set(tumour_gt) - set(reference_nt)
        variant_nt = variant_set.pop() if len(variant_set)!=0 else ''
        
        if variant_nt=='':
            print('Warning: no valid genotypes for variant %s:%d; skipping.'%(record.CHROM,record.POS))
            continue
            
        normal = record.genotype('NORMAL')
        tumour = record.genotype('TUMOUR')

        tumor_reads = {
          'forward': {
            'A': int(tumour['FAZ']),
            'C': int(tumour['FCZ']),
            'G': int(tumour['FGZ']),
            'T': int(tumour['FTZ']),
          },
          'reverse': {
            'A': int(tumour['RAZ']),
            'C': int(tumour['RCZ']),
            'G': int(tumour['RGZ']),
            'T': int(tumour['RTZ']),
          },
        }

        ref_reads = tumor_reads['forward'][reference_nt] + tumor_reads['reverse'][reference_nt]
        variant_reads = tumor_reads['forward'][variant_nt] + tumor_reads['reverse'][variant_nt]
        total_reads = ref_reads + variant_reads
        
        if variant_reads!=0:
            tmp = np.array((record.CHROM,record.POS,'',ref_reads,variant_reads),dtype=snv_dtype)
            snv_df = np.append(snv_df,tmp)

    #import matplotlib.pyplot as plt
    #fig, axes = plt.subplots(1, 1, sharex=False, sharey=False)
    #dep = snv_df['ref']+snv_df['var']
    #sup = map(float,snv_df['var'])
    #axes.hist(sup/dep);plt.savefig('/home/mcmero/Desktop/test')
    
    return pd.DataFrame(snv_df)
