'''
Run clustering and tree building on sample inputs
'''
import os
import numpy as np
import ipdb
import re
import pandas as pd
import vcf
from operator import methodcaller
from . import run_clus
from . import parameters as params

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

def run_cnv_filter(df_flt,cnv,neutral=False,perc=85):
    '''
    filter based on either CNV neutral, or presence of CNV vals
    '''
    if len(cnv)>0 and neutral:
        # filter out copy-aberrant SVs and outying norm read counts (>1-percentile)
        # major and minor copy-numbers must be 1
        gt1_is_neutral = map(is_clonal_neutral,df_flt.gtype1.values)
        gt2_is_neutral = map(is_clonal_neutral,df_flt.gtype2.values)
        is_neutral = np.logical_and(gt1_is_neutral,gt2_is_neutral)
        df_flt = df_flt[is_neutral & (df_flt.norm_mean<np.percentile(df.norm_mean,perc))] 
    elif len(cnv)>0:
        #TODO: would be good to make a filter to check if normal reads differ wildly from possible CNV states
        df_flt = df_flt[np.logical_or(df_flt.gtype1.values!='',df_flt.gtype2.values!='')]
    return df_flt

def cnv_overlaps(bp_pos,cnv_start,cnv_end):
    return (bp_pos >= cnv_start and bp_pos <= cnv_end)

#def find_cn_cols(bp_chr,bp_pos,cnv,cols=['nMaj1_A','nMin1_A','frac1_A','nMaj2_A','nMin2_A','frac2_A']):
#    for idx,c in cnv.iterrows():
#        if str(bp_chr)==c['chr'] and cnv_overlaps(bp_pos,c['startpos'],c['endpos']):
#            return c[cols[0]],c[cols[1]],c[cols[2]],c[cols[3]],c[cols[4]],c[cols[5]]
#    return float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')

def find_cn_col(bp_chr,bp_pos,cnv):
    for idx,c in cnv.iterrows():
        if str(bp_chr)==c['chr'] and cnv_overlaps(bp_pos,c['startpos'],c['endpos']):
            return c
    return pd.Series()

def get_genotype(cnrow):
    gtype = []
    for col in cnrow.iteritems():
        mj = re.search("n(Maj)(?P<fraction>\d_\w)",col[0])
        if len(gtype)==2: break
        if mj:
            if np.isnan(col[1]) or col[1]==0: break
            nmin = 'nMin' + mj.group("fraction") #corresponding minor allele
            frac = 'frac' + mj.group("fraction") #allele fraction
            gtype.append([col[1],cnrow[nmin],cnrow[frac]])
    g_str = ["%d,%d,%f"%(gt[0],gt[1],gt[2]) for gt in gtype]
    return '|'.join(g_str)

def load_snvs(snvs):
    vcf_reader = vcf.Reader(filename=snvs)
    snv_dtype = [('chrom','S50'),('pos',int),('gtype','S50'),('ref',int),('var',int)]
    snv_df = np.empty([0,5],dtype=snv_dtype)
    sample = vcf_reader.samples[0]

    if len(vcf_reader.samples)>1 or sample=='':
        print('Skipping SNV clustering, VCF does not have exactly one sample')
        return pd.DataFrame()
        
    for record in vcf_reader:
        ad = record.genotype(sample)['AD']
        if ad[1]!=0:
            tmp = np.array((record.CHROM,record.POS,'',ad[0],ad[1]),dtype=snv_dtype)
            snv_df = np.append(snv_df,tmp)

    return pd.DataFrame(snv_df)

def is_same_sv(sv1,sv2):
    sv1_chr1, sv1_bp1, sv1_chr2, sv1_bp2 = sv1
    sv2_chr1, sv2_bp1, sv2_chr2, sv2_bp2 = sv2

    if sv1_chr1==sv2_chr1 and sv1_chr2==sv2_chr2:
        if abs(sv1_bp1-sv2_bp1)<params.gl_th and abs(sv1_bp2-sv2_bp2)<params.gl_th:
            return True
    return False

def run(samples,svs,gml,cnvs,rlens,inserts,pis,ploidies,out,n_runs,num_iters,burn,thin,beta,neutral,snvs):
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
            snv_df = load_snvs(snvs)

        filt_svs_file = '%s/%s_filtered_svs.tsv'%(out,sample)
        if os.path.isfile(filt_svs_file):
            print('Found filtered SVs file, running clustering on these variants')
            sv_df = pd.read_csv(filt_svs_file,delimiter='\t',na_values='')
            sv_df = pd.DataFrame(sv_df).fillna('')
            print('Clustering with %d SVs and %d SNVs'%(len(sv_df),len(snv_df)))
            run_clus.infer_subclones(sample,sv_df,pi,rlen,insert,ploidy,out,n_runs,num_iters,burn,thin,beta,snv_df)
            return None
        else:
            dat = pd.read_csv(sv,delimiter='\t')
            df = pd.DataFrame(dat)
            df['norm_mean'] = map(np.mean,zip(df['norm1'].values,df['norm2'].values))
            sv_df = run_simple_filter(df,rlen,insert)

        if gml!="":
            df_gml = pd.DataFrame(pd.read_csv(gml,delimiter='\t'))
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
            sv_df = sv_df.drop(germline,axis=0)

        if cnv!="":
            cnv = pd.read_csv(cnv,delimiter='\t')
            cnv_sv_df = pd.DataFrame(cnv)
            cnv_sv_df['chr'] = map(str,cnv_sv_df['chr'])
           
            if len(snv_df)>0:
                for idx,snv in snv_sv_df.iterrows():
                    cn = find_cn_col(snv.chrom,snv.pos,cnv)
                    snv_sv_df['gtype'][idx] = get_genotype(cn)
                snv_sv_df = snv_sv_df[snv_sv_df.gtype!='']

            sv_df['gtype1'] = "" 
            sv_df['gtype2'] = "" 
            for idx,d in sv_df.iterrows():
                cn1 = find_cn_col(d.bp1_chr,d.bp1_pos,cnv)
                cn2 = find_cn_col(d.bp2_chr,d.bp2_pos,cnv)
                sv_df['gtype1'][idx] = get_genotype(cn1)
                sv_df['gtype2'][idx] = get_genotype(cn2)

            sv_df = run_cnv_filter(sv_df,cnv,neutral)
        else:
            print('No Battenberg input defined, assuming all loci are copy-number neutral')
            sv_df['gtype1'] = '1,1,1.000000'
            sv_df['gtype2'] = '1,1,1.000000'
            sv_df = sv_df[(sv_df.norm_mean<np.percentile(sv_df.norm_mean,85))]                

            
        sv_df.index = range(len(sv_df))
        sv_df.to_csv('%s/%s_filtered_svs.tsv'%(out,sample),sep='\t',index=False,na_rep='')
        
        with open('%s/purity_ploidy.txt'%out,'w') as outf:
            outf.write("sample\tpurity\tploidy\n")
            outf.write('%s\t%f\t%f\n'%(sample,pi,ploidy))
      
        print('Clustering with %d SVs and %d SNVs'%(len(sv_df),len(snv_df)))
        run_clus.infer_subclones(sample,sv_df,pi,rlen,insert,ploidy,out,n_runs,num_iters,burn,thin,beta,snv_df)

