'''
Run clustering and tree building on sample inputs
'''

import numpy as np
import phylo_sv as ps
import ipdb
import pandas as pd

def run_filter(df,rlen,insert,cnv=np.empty(0),perc=85):
    #filter based on fragment length
    itx = df['bp1_chr']!=df['bp2_chr']
    frag_len = 2*rlen+insert
    df_flt = df[ itx | (abs(df['bp1_pos']-df['bp2_pos'])>frag_len) ]
    df_flt = df_flt[df_flt['support']>1]
    df['norm_mean'] = map(np.mean,zip(df['norm1'].values,df['norm2'].values))
    if len(cnv)>0:
        df_flt = df[(df.bp1_frac1A==1) & (df.bp2_frac1A==1) & (df.bp1_maj_cnv+df.bp1_min_cnv<=3) & \
            (df.bp2_maj_cnv+df.bp2_min_cnv<=3) & (df.norm_mean<np.percentile(df.norm_mean,perc))]
    return df_flt

def cnv_overlaps(bp_pos,cnv_start,cnv_end):
    return (bp_pos >= cnv_start and bp_pos <= cnv_end)

def find_cn_cols(bp_chr,bp_pos,cnv,cols=['nMaj1_A','nMin1_A','frac1_A']):
    for idx,c in cnv.iterrows():
        if str(bp_chr)==c['chr'] and cnv_overlaps(bp_pos,c['startpos'],c['endpos']):
            return c[cols[0]],c[cols[1]],c[cols[2]]
    return float('nan'),float('nan'),float('nan')

def run(samples,svs,gml,cnvs,rlens,inserts,pis,out):
    if gml!="":
        #run germline filtering
        print
   
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
           
            # set BattenBerg fields
            df['bp1_maj_cnv'], df['bp1_min_cnv'], df['bp1_frac1A'] = float('nan'), float('nan'), float('nan')
            df['bp2_maj_cnv'], df['bp2_min_cnv'], df['bp2_frac1A'] = float('nan'), float('nan'), float('nan')

            for idx,d in df.iterrows():
                df['bp1_maj_cnv'][idx], df['bp1_min_cnv'][idx], df['bp1_frac1A'][idx] = \
                        find_cn_cols(d['bp1_chr'],d['bp1_pos'],cnv)
                df['bp2_maj_cnv'][idx], df['bp2_min_cnv'][idx], df['bp2_frac1A'][idx] = \
                        find_cn_cols(d['bp2_chr'],d['bp2_pos'],cnv)
        
        df_flt = run_filter(df,rlen,insert,cnv)
        ps.infer_tree(df_flt,pi,rlen,insert,out)
