'''
Run clustering and tree building on sample inputs
'''
import os
import numpy as np
import ipdb
import re
import pandas as pd
import vcf
# import cProfile, pstats, StringIO 
from operator import methodcaller
from . import run_clus
from . import cluster
from . import parameters as params

def get_outlier_ranges(vals):
    q1 = np.percentile(vals,25)
    q3 = np.percentile(vals,75)
    iqr = max(q3 - q1, 1)
    lower = q1 - 1.5*iqr
    upper = q3 + 1.5*iqr
    return [lower,upper]

def run_simple_filter(df,rlen,insert):
    '''
    filter based on presence of supporting reads and > fragment length
    '''
    span = np.array(df.spanning)
    split = np.array(df.bp1_split+df.bp2_split)
    df_flt = df[(span+split)>=1]

    itx = df_flt['bp1_chr']!=df_flt['bp2_chr']
    frag_len = 2*rlen+insert
    df_flt = df_flt[ itx | (abs(df_flt['bp1_pos']-df_flt['bp2_pos'])>2*rlen) ]
    
    return df_flt

def is_clonal_neutral(gtype):
    if gtype=='': return False
    gtype = map(methodcaller('split',','),gtype.split('|'))
    if len(gtype)==1:
        gt = map(float,gtype[0])
        return (gt[0]==1 and gt[1]==1 and gt[2]==1)
    return False

def filter_outlying_norm_wins(df_flt):
    bp1_win, bp2_win = cluster.normalise_wins_by_cn(df_flt)
    bp1_ranges = get_outlier_ranges(bp1_win)
    bp2_ranges = get_outlier_ranges(bp2_win)

    bp1_flt = np.logical_and(bp1_win>bp1_ranges[0],bp1_win<bp1_ranges[1])
    bp2_flt = np.logical_and(bp2_win>bp2_ranges[0],bp2_win<bp2_ranges[1])
    df_flt = df_flt[np.logical_and(bp1_flt,bp2_flt)]
    
    return df_flt

def run_cnv_filter(df_flt,cnv,neutral=False,are_snvs=False):
    '''
    filter based on either CNV neutral, or presence of CNV vals
    '''
    if len(cnv)>0 and neutral:
        # filter out copy-aberrant SVs and outying norm read counts (>1-percentile)
        # major and minor copy-numbers must be 1
        is_neutral = []
        if are_snvs:
            is_neutral = map(is_clonal_neutral,df_flt.gtype.values)
            df_flt = df_flt[is_neutral]
            depths = df_flt['ref'].values+df_flt['var'].values
            dep_ranges = get_outlier_ranges(depths)
            df_flt = df_flt[np.logical_and(depths>dep_ranges[0],depths<dep_ranges[1])]
        else: 
            gt1_is_neutral = map(is_clonal_neutral,df_flt.gtype1.values)
            gt2_is_neutral = map(is_clonal_neutral,df_flt.gtype2.values)
            is_neutral = np.logical_and(gt1_is_neutral,gt2_is_neutral)
            df_flt = df_flt[is_neutral] 
            #df_flt = df_flt[is_neutral & (df_flt.norm_mean<np.percentile(df_flt.norm_mean,perc))] 
            df_flt = filter_outlying_norm_wins(df_flt)
    elif len(cnv)>0:
        if are_snvs:
            df_flt = df_flt[df_flt.gtype.values!='']

            # weight ranges by copy-numbers
            depths = df_flt['ref'].values + df_flt['var'].values
            cns = cluster.get_weighted_cns(df_flt.gtype.values)
            cn_nonzero = np.logical_not(cns==0) 
            depths[cn_nonzero] = (depths/cns)[cn_nonzero]

            dep_ranges = get_outlier_ranges(depths)
            df_flt = df_flt[np.logical_and(depths>dep_ranges[0],depths<dep_ranges[1])]
        else:
            df_flt = df_flt[np.logical_or(df_flt.gtype1.values!='',df_flt.gtype2.values!='')]
            df_flt = filter_outlying_norm_wins(df_flt)
            #df_flt = df_flt[np.logical_and(df_flt.bp1_win_norm<600,df_flt.bp2_win_norm<600)]
    if are_snvs:
        df_flt = df_flt.fillna('')
        df_flt = df_flt[np.logical_and(df_flt.gtype.values!='',df_flt.chrom.values!='')]
    return df_flt

#def cnv_overlaps(bp_pos,cnv_start,cnv_end):
#    return (bp_pos >= cnv_start and bp_pos <= cnv_end)

#def find_cn_cols(bp_chr,bp_pos,cnv,cols=['nMaj1_A','nMin1_A','frac1_A','nMaj2_A','nMin2_A','frac2_A']):
#    for idx,c in cnv.iterrows():
#        if str(bp_chr)==c['chr'] and cnv_overlaps(bp_pos,c['startpos'],c['endpos']):
#            return c[cols[0]],c[cols[1]],c[cols[2]],c[cols[3]],c[cols[4]],c[cols[5]]
#    return float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')

#def find_cn_col(bp_chr,bp_pos,cnv):
#    for idx,c in cnv.iterrows():
#        if str(bp_chr)==c['chr'] and cnv_overlaps(bp_pos,c['startpos'],c['endpos']):
#            return c
#    return pd.Series()
#
#def get_genotype(cnrow):
#    #print('get genotype for %s:%d' % (cnrow['chr'],cnrow['startpos'])) 
#    gtype = []
#    bb_fields = ['nMaj1_A','nMin1_A','frac1_A','nMaj2_A','nMin2_A','frac2_A']
#    maj1A, min1A, frac1A = cnrow[bb_fields[:3]].values 
#    
#    if np.nan in [maj1A, min1A, frac1A]:
#        raise('Standard battenberg fields not found. Please check your copy-number input.')
#
#    g_str = '%d,%d,%f' % (maj1A, min1A, frac1A)
#    if frac1A!=1:
#        maj2A, min2A, frac2A = cnrow[bb_fields[3:]].values
#        g_str += '|%d,%d,%f' % (maj2A, min2A, frac2A)
#
#    return g_str

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

def match_copy_numbers(var_df, cnv_df, bp_fields=['bp1_chr','bp1_pos'], gtype_field='gtype1'):
    
#    var_df['gtype1'] = ''
#    var_df['gtype2'] = ''
#    sv_arr = np.array(zip(var_df.index.values,
#                          var_df['bp1_chr'].values,
#                          var_df['bp1_pos'].values,
#                          var_df['bp2_chr'].values,
#                          var_df['bp2_pos'].values),
#                     dtype=[('idx',int),('bp1_chr','S50'),\
#                            ('bp1_pos',int),('bp2_chr','S50'),('bp2_pos',int)])
#
#    for sv in sv_arr:
#        idx = sv['idx']
#        for cnv in cnv_arr:
#            if sv['bp1_chr'] == cnv['chr'] and cnv_overlaps(sv['bp1_pos'],cnv['startpos'],cnv['endpos']):
#               var_df.loc[idx,'gtype1'] = cnv['gtype']
#            if sv['bp2_chr'] == cnv['chr'] and cnv_overlaps(sv['bp2_pos'],cnv['startpos'],cnv['endpos']):
#               var_df.loc[idx,'gtype2'] = cnv['gtype']
#
#    return var_df 
    var_df[gtype_field] = ''
    chrom_field, pos_field = bp_fields
    var_df[chrom_field] = map(str,var_df[chrom_field].values)
    var_df[pos_field] = map(int,var_df[pos_field].values)
    bp_chroms = np.unique(var_df[chrom_field].values)
    bp_chroms = sorted(bp_chroms, key=lambda item: (int(item.partition(' ')[0]) \
                        if item[0].isdigit() else float('inf'), item))
    
    for bchr in bp_chroms:
        print('Matching copy-numbers for chrom %s'%bchr)
        gtypes = []
        current_chr = var_df[chrom_field].values==bchr
        sv_tmp = var_df[current_chr]
        cnv_tmp = cnv_df[cnv_df['chr']==bchr]
        
        if len(cnv_tmp)==0:
            continue

        for pos in sv_tmp[pos_field].values:
            #print("Checking overlaps for pos %d\n" %pos)
            overlaps = np.logical_and(pos >= cnv_tmp.startpos.values, pos <= cnv_tmp.endpos.values)
            match = cnv_tmp[overlaps]
            gtype = match.loc[match.index[0]].gtype if len(match)>0 else '' 
            gtypes.append(gtype)
        
        var_indexes = var_df[current_chr].index.values
        var_df.loc[var_indexes,gtype_field] = gtypes
    return var_df

def load_svs(sv_file, rlen, insert):
    dat = pd.read_csv(sv_file,delimiter='\t',dtype=None, low_memory=False)
    sv_df = pd.DataFrame(dat)
    sv_df['norm_mean'] = map(np.mean,zip(sv_df['norm1'].values,sv_df['norm2'].values))
    return run_simple_filter(sv_df,rlen,insert)

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
    #TODO: remove sample
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
    #ipdb.set_trace()
    return pd.DataFrame(snv_df)

def is_same_sv(sv1,sv2):
    sv1_chr1, sv1_bp1, sv1_chr2, sv1_bp2 = sv1
    sv2_chr1, sv2_bp1, sv2_chr2, sv2_bp2 = sv2

    if sv1_chr1==sv2_chr1 and sv1_chr2==sv2_chr2:
        if abs(sv1_bp1-sv2_bp1)<params.gl_th and abs(sv1_bp2-sv2_bp2)<params.gl_th:
            return True
    return False

def filter_germline(gml_file,sv_df,rlen,insert):
    #TODO: optimise (avoid using iterrows)
    print("Filtering out germline SVs...")
    df_gml = pd.DataFrame(pd.read_csv(gml_file,delimiter='\t',dtype=None,low_memory=False))
    df_gml = run_simple_filter(df_gml,rlen,insert)
    germline = []
    
    for idx_sv,sv in sv_df.iterrows():
        sv_loc = [str(sv.bp1_chr),int(sv.bp1_pos),str(sv.bp2_chr),int(sv.bp2_pos)]
        same_chrs = np.logical_and(df_gml.bp1_chr.values==sv_loc[0],df_gml.bp2_chr.values==sv_loc[2])
        df_gml_tmp = df_gml[same_chrs]
        for idx_gml,sv_gml in df_gml_tmp.iterrows():
            sv_loc_gml = [str(sv_gml.bp1_chr),int(sv_gml.bp1_pos),str(sv_gml.bp2_chr),int(sv_gml.bp2_pos)]
            if sv_gml.support>0 and is_same_sv(sv_loc,sv_loc_gml):
                germline.append(idx_sv)
                break
                
    print("Found %d SVs in the germline!" % len(germline))
    return sv_df.drop(germline,axis=0)

def run(samples,svs,gml,cnvs,rlens,inserts,pis,ploidies,out,n_runs,num_iters,burn,thin,beta,neutral,snvs,snv_format,merge_clusts,use_map):

    # pr = cProfile.Profile()
    # pr.enable()
    
    if not os.path.exists(out):
        os.makedirs(out)

    if len(samples)>1:
        print("Multiple sample processing not yet implemented")
        #TODO: processing of multiple samples
        #for sm,sv,cnv,rlen,insert,pi in zip(samples,svs,cnvs,rlens,inserts,pis):
        #    dat = pd.read_csv(sv,delimiter='\t')
        #    df = pd.DataFrame(dat)
        #    svinfo = pd.DataFrame(dat)
    else:
        sample,sv,cnv,rlen,insert,pi,ploidy = samples[0],svs[0],cnvs[0],rlens[0],inserts[0],pis[0],ploidies[0]
        sv_df  = pd.DataFrame()
        cnv_df = pd.DataFrame()
        snv_df = pd.DataFrame()
       
        if snvs!="":
            if snv_format == 'sanger':
                snv_df = load_snvs_sanger(snvs)
            elif snv_format == 'mutect':
                snv_df = load_snvs_mutect(snvs,sample)
            elif snv_format == 'mutect_callstats':
                snv_df = load_snvs_mutect_callstats(snvs)

        filt_svs_file = '%s/%s_filtered_svs.tsv'%(out,sample)
        sv_df = load_svs(sv,rlen,insert)
#        if os.path.isfile(filt_svs_file):
#            print('Found filtered SVs file, running clustering on these variants')
#            sv_df = pd.read_csv(filt_svs_file,delimiter='\t',na_values='')
#            sv_df = pd.DataFrame(sv_df).fillna('')
#            print('Clustering with %d SVs and %d SNVs'%(len(sv_df),len(snv_df)))
#            run_clus.infer_subclones(sample,sv_df,pi,rlen,insert,ploidy,out,n_runs,num_iters,burn,thin,beta,snv_df)
#            return None

        if gml!="":
            sv_df = filter_germline(gml,sv_df,rlen,insert)
        
        if cnv!="":
            cnv_df = load_cnvs(cnv)
            sv_df = match_copy_numbers(sv_df,cnv_df) 
            sv_df = match_copy_numbers(sv_df,cnv_df,['bp2_chr','bp2_pos'],'gtype2') 
            sv_df = run_cnv_filter(sv_df,cnv,neutral)
            
            if len(snv_df)>0:
                snv_df = match_copy_numbers(snv_df,cnv_df,['chrom','pos'],'gtype')
                snv_df = run_cnv_filter(snv_df,cnv,neutral,are_snvs=True)
        else:
            print('No Battenberg input defined, assuming all loci are copy-number neutral')
            sv_df['gtype1'] = '1,1,1.000000'
            sv_df['gtype2'] = '1,1,1.000000'
            #sv_df = sv_df[(sv_df.norm_mean<np.percentile(sv_df.norm_mean,85))]                
            sv_df = filter_outlying_norm_wins(sv_df)

        #reindex data frames
        sv_df.index = range(len(sv_df))
        if len(snv_df)>0: snv_df.index = range(len(snv_df))

        sv_df.index = range(len(sv_df))
        sv_df.to_csv('%s/%s_filtered_svs.tsv'%(out,sample),sep='\t',index=False,na_rep='')
        
        with open('%s/purity_ploidy.txt'%out,'w') as outf:
            outf.write("sample\tpurity\tploidy\n")
            outf.write('%s\t%f\t%f\n'%(sample,pi,ploidy))
      
        if len(sv_df) < 5:
            print("Less than 5 post-filtered SVs. Clustering not recommended for this sample. Exiting.")
        else:
            print('Clustering with %d SVs and %d SNVs'%(len(sv_df),len(snv_df)))        
            run_clus.infer_subclones(sample,sv_df,pi,rlen,insert,ploidy,out,n_runs,num_iters,burn,thin,beta,snv_df,merge_clusts,use_map)

        # pr.disable()
        # s = StringIO.StringIO()
        # sortby = 'cumulative'
        # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        # ps.print_stats()
        # print s.getvalue()

