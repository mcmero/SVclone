'''
Run clustering and tree building on sample inputs
'''
import os
import sys
import numpy as np
import ipdb
import pandas as pd
from . import build_phyl

def run_filter(df,rlen,insert,cnv,ploidy,perc=85):
    df_flt = df[df['support']>1]
    
    #filter based on fragment length
    itx = df_flt['bp1_chr']!=df_flt['bp2_chr']
    frag_len = 2*rlen+insert
    df_flt = df_flt[ itx | (abs(df_flt['bp1_pos']-df_flt['bp2_pos'])>frag_len) ]
    
    if len(cnv)>0:
        # filter out copy-aberrant SVs and outying norm read counts (>1-percentile)
        #cn1 = np.array(map(round,df.bp1_maj_cnv+df.bp1_min_cnv))
        #cn2 = np.array(map(round,df.bp2_maj_cnv+df.bp2_min_cnv))
        #df_flt = df[(df.bp1_frac1A==1) & (df.bp2_frac1A==1) & (cn1==2) & (cn2==2) & \
        #        (df.norm_mean<np.percentile(df.norm_mean,perc))]
        # major and minor copy-numbers must be 1
        cn1_maj = np.array(map(round,df.bp1_maj_cnv))
        cn1_min = np.array(map(round,df.bp1_min_cnv))
        cn2_maj = np.array(map(round,df.bp2_maj_cnv))
        cn2_min = np.array(map(round,df.bp2_min_cnv))
        bp1_cnn = np.logical_and(cn1_maj==1,cn1_min==1)
        bp2_cnn = np.logical_and(cn2_maj==1,cn2_min==1)
        df_flt = df[(df.bp1_frac1A==1) & bp1_cnn & bp2_cnn & \
                    (df.norm_mean<np.percentile(df.norm_mean,perc))]
    df_flt.index = range(len(df_flt))
    return df_flt

def cnv_overlaps(bp_pos,cnv_start,cnv_end):
    return (bp_pos >= cnv_start and bp_pos <= cnv_end)

def find_cn_cols(bp_chr,bp_pos,cnv,cols=['nMaj1_A','nMin1_A','frac1_A']):
    for idx,c in cnv.iterrows():
        if str(bp_chr)==c['chr'] and cnv_overlaps(bp_pos,c['startpos'],c['endpos']):
            return c[cols[0]],c[cols[1]],c[cols[2]]
    return float('nan'),float('nan'),float('nan')

def run(samples,svs,gml,cnvs,rlens,inserts,pis,ploidy,out,num_iters):
    if not os.path.exists(out):
        os.makedirs(out)

    if gml!="":
        #TODO: implement germline filtering
        print("Germline filtering not yet implemented")
   
    if len(samples)>1:
        print("Multiple sample processing not yet implemented")
        #TODO: processing of multiple samples
        #for sm,sv,cnv,rlen,insert,pi in zip(samples,svs,cnvs,rlens,inserts,pis):
        #    dat = pd.read_csv(sv,delimiter='\t')
        #    df = pd.DataFrame(dat)
        #    svinfo = pd.DataFrame(dat)
    else:
        sv,cnv,rlen,insert,pi = svs[0],cnvs[0],rlens[0],inserts[0],pis[0]
        dat = pd.read_csv(sv,delimiter='\t')
        df = pd.DataFrame(dat)
        cnv_df = pd.DataFrame()

        if cnv!="":
            cnv = pd.read_csv(cnv,delimiter='\t')
            cnv_df = pd.DataFrame(cnv)
            cnv_df['chr'] = map(str,map(int,cnv_df['chr']))
           
            # set Battenberg fields
            df['bp1_maj_cnv'], df['bp1_min_cnv'], df['bp1_frac1A'] = float('nan'), float('nan'), float('nan')
            df['bp2_maj_cnv'], df['bp2_min_cnv'], df['bp2_frac1A'] = float('nan'), float('nan'), float('nan')

            for idx,d in df.iterrows():
                df['bp1_maj_cnv'][idx], df['bp1_min_cnv'][idx], df['bp1_frac1A'][idx] = \
                        find_cn_cols(d['bp1_chr'],d['bp1_pos'],cnv)
                df['bp2_maj_cnv'][idx], df['bp2_min_cnv'][idx], df['bp2_frac1A'][idx] = \
                        find_cn_cols(d['bp2_chr'],d['bp2_pos'],cnv)
            df = df[pd.notnull(df['bp1_frac1A'])]        

        df['norm_mean'] = map(np.mean,zip(df['norm1'].values,df['norm2'].values))
        df_flt = run_filter(df,rlen,insert,cnv,ploidy)
        if len(df_flt) < 5:
            print("Less than 5 post-filtered SVs. Clustering not recommended for this sample. Exiting.")
            return None
        build_phyl.infer_subclones(df_flt,pi,rlen,insert,out,num_iters)
