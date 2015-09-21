import pandas as pd
import numpy as np
import vcf

from . import cluster

def get_sv_vals(sv_df,no_adjust):
    combos = sv_df.apply(cluster.get_sv_allele_combos,axis=1)
    sides = sv_df.preferred_side.values
    cn_states = [cn[side] for cn,side in zip(combos,sides)]
    if no_adjust:
        sup = sv_df.adjusted_support.map(float).values
        dep = sv_df.adjusted_depth.map(float).values
        Nvar = len(sv_df)
        return sup,dep,cn_states,Nvar
    else:
        sup  = sv_df.support.map(float).values
        norm = zip(sv_df.norm1.values,sv_df.norm2.values)
        norm = np.array([float(n[side]) for n,side in zip(norm,sv_df.preferred_side.values)])
        dep  = norm+sup 
        Nvar = len(sv_df)
        return sup,dep,cn_states,Nvar

def get_snv_vals(df):
    n = df['ref'].map(float).values
    b = df['var'].map(float).values

    def get_snv_allele_combos(snv):
        return cluster.get_allele_combos(snv.gtype.split('|'))
     
    combos = df.apply(get_snv_allele_combos,axis=1)
    return b,(n+b),combos,len(b)

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

#def get_sv_vals(df,rlen,pi,pl):
#    n = zip(np.array(df.norm1.values),np.array(df.norm2.values))
#    s = np.array(df.bp1_split.values+df.bp2_split.values)
#    d = np.array(df.spanning.values)
#    #cn_r,cn_v,mu_v = [],[],[]
#    
#    combos = df.apply(get_sv_allele_combos,axis=1)
#    #cn_r,cn_v,mu_v = combos[0],combos[1],combos[2]
#
##    for idx,sv in df.iterrows():
##        cn_tmp = tuple([tuple(sv.gtype1.split('|')),tuple(sv.gtype2.split('|'))])
##        cnr_bp1,cnv_bp1,mu_bp1 = get_allele_combos_tuple(cn_tmp[0])
##        cnr_bp2,cnv_bp2,mu_bp2 = get_allele_combos_tuple(cn_tmp[1])
##        cn_r.append(tuple([cnr_bp1,cnr_bp2]))
##        cn_v.append(tuple([cnv_bp1,cnv_bp2]))
##        mu_v.append(tuple([mu_bp1,mu_bp2]))
#
#    sup = np.array(d+s,dtype=float)
#    n_bp1,n_bp2 = [ni[0] for ni in n],[ni[1] for ni in n]
#    Nvar = len(sup)
#    
#    # calculate the average depth using windowed counts
#    #win1 = (df.bp1_win_norm.values*rlen)/(param.window*2)
#    #win2 = (df.bp2_win_norm.values*rlen)/(param.window*2)
#    #av_cov = np.mean([win1,win2])
#     
#    sides = np.zeros(Nvar,dtype=int)
#    sides[df.gtype1.values=='']=1 #no genotype for bp1, use bp2
#
#    #has_both_gts = np.logical_and(df.gtype1.values!='',df.gtype2.values!='')
#    
#    #bp1_win, bp2_win = normalise_wins_by_cn(df)
#    #bp1_win, bp2_win = (bp2_win*rlen)/(param.window*2), (bp2_win*rlen)/(param.window*2)
#
#    # for sides with gtypes for both sides, pick the side where the adjusted window count is closest to the average coverage
#    #av_cov = np.mean([bp1_win,bp2_win])
#    #dev_from_cov = np.array(zip(abs(bp1_win-av_cov),abs(bp2_win-av_cov)))
#
#    # ALT: pick the side where the norm count is closest to the mean coverage
#    #dev_from_cov = np.array(zip(abs(n_bp1-av_cov),abs(n_bp2-av_cov)))
#    
#    # prefer sides with subclonal genotype data    
#    gt1_sc = np.array(map(len,map(methodcaller("split","|"),df.gtype1.values)))>1
#    gt2_sc = np.array(map(len,map(methodcaller("split","|"),df.gtype1.values)))>1    
#    one_sc = np.logical_xor(gt1_sc,gt2_sc)
#
#    exclusive_subclones = zip(df.gtype1.values[one_sc],df.gtype2.values[one_sc]) 
#    sides[one_sc] = [0 if gt1!='' else 1 for gt1,gt2 in exclusive_subclones]
#
#    has_both_gts = np.logical_and(df.gtype1.values!='',df.gtype2.values!='')
#    #sides[has_both_gts] = 
#    #sides[has_both_gts] = [np.where(x==min(x))[0][0] for x in dev_from_cov[has_both_gts]]
#
#    norm = np.array([ni[si] for ni,si in zip(n,sides)])
#    #norm = map(np.mean,n)
#
#    # both sides have genotypes, both either subclonal or clonal
#    # in this case, just take the simple normal mean of the two sides
#    both_gts_same_type = np.logical_and(has_both_gts,one_sc==False)
#    norm[both_gts_same_type] = map(np.mean,np.array(n)[both_gts_same_type])
#
#    #return sup,dep,cn_r,cn_v,mu_v,sides,Nvar
#    cn_states = [cn[side] for cn,side in zip(combos,sides)]
#
#    dep = np.array(norm+sup,dtype=float)
#    return sup,dep,cn_states,Nvar
