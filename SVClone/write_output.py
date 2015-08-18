import numpy as np
import pandas as pd
import os
import ipdb
from operator import methodcaller

from . import cluster

def dump_trace(clus_info,center_trace,outf):
    traces = np.array([center_trace[:,cid] for cid in clus_info.clus_id.values])
    df_traces = pd.DataFrame(np.transpose(traces),columns=clus_info.clus_id)
    df_traces.to_csv(outf,sep='\t',index=False)

def write_out_files_snv(df,clus_info,clus_members,df_probs,clus_cert,clus_out_dir,sample,pi):
    clus_out_dir = '%s/snvs'%clus_out_dir
    if not os.path.exists(clus_out_dir):
        os.makedirs(clus_out_dir)

    clus_info.to_csv('%s/clusters.txt'%(clus_out_dir),sep='\t',index=False)
    with open('%s/number_of_clusters.txt'%clus_out_dir,'w') as outf:
        outf.write("sample\tclusters\n")
        outf.write('%s\t%d\n'%(sample,len(clus_info)))
    
    cn_dtype =  [('chr1','S50'),
                 ('pos1',int),
                 ('total_copynumber1',int),
                 ('no_chrs_bearing_mutation2',int)]

    mlcn_dtype =[('chrom','S50'),
                 ('pos',int),
                 ('bb_CN','S50'),
                 ('most_likely_ref_copynumber',int),
                 ('most_likely_variant_copynumber',int),
                 ('prop_chrs_bearing_mutation',float)]

    cmem = np.hstack(clus_members)
    cn_vect = np.empty((0,len(cmem)),dtype=cn_dtype)
    mlcn_vect = np.empty((0,len(cmem)),dtype=mlcn_dtype)
    #clus_snvs = df.loc[cmem].copy()
    df = df.fillna('')
    df = df[df.chrom.values!='']

#    sup,n,cn_r,cn_v,mu_v = cluster.get_snv_vals(df)
#    dep = sup + n
#    phis = clus_cert.average_ccf.values    
#
#    for idx,snv in df.iterrows():
#        gtype = snv['gtype'].split('|')
#        gtype = map(methodcaller('split', ','), gtype)
#
#        maj_cn,min_cn = map(float,gtype[0])[:2] if gtype[0][0]!='' else [0.,0.]
#        tot_cn = maj_cn+min_cn
#
#        chrom = str(snv['chrom'])
#        pos = int(snv['pos'])
#        ref_cn, sc_cn, freq = cluster.get_most_likely_cn(cn_r[idx][0],cn_v[idx][0],\
#                                                         mu_v[idx][0],sup[idx],dep[idx],phis[idx],pi) 
#        
#        cn_new_row = np.array([(chrom,pos,tot_cn,int(sc_cn*freq))],dtype=cn_dtype)
#        cn_vect = np.append(cn_vect,cn_new_row)
#        
#        ml_new_row = np.array([(chrom,pos,snv['gtype'],ref_cn,sc_cn,freq)],dtype=mlcn_dtype)
#        mlcn_vect = np.append(mlcn_vect,ml_new_row)
#        
#    pd.DataFrame(mlcn_vect).to_csv('%s/%s_most_likely_copynumbers.txt'%(clus_out_dir,sample),sep='\t',index=False)
#    pd.DataFrame(cn_vect).to_csv('%s/%s_copynumber.txt'%(clus_out_dir,sample),sep='\t',index=False)
#    df_probs.to_csv('%s/%s_assignment_probability_table.txt'%(clus_out_dir,sample),sep='\t',index=False)
#    clus_cert.to_csv('%s/%s_cluster_certainty.txt'%(clus_out_dir,sample),sep='\t',index=False)

def write_out_files(df,clus_info,clus_members,df_probs,clus_cert,clus_out_dir,sample,pi,rlen):
    
    clus_info.to_csv('%s/clusters.txt'%(clus_out_dir),sep='\t',index=False)
    with open('%s/number_of_clusters.txt'%clus_out_dir,'w') as outf:
        outf.write("sample\tclusters\n")
        outf.write('%s\t%d\n'%(sample,len(clus_info)))
    
    cn_dtype =  [('chr1','S50'),
                 ('pos1',int),
                 ('total_copynumber1',int),
                 ('no_chrs_bearing_mutation1',int),
                 ('chr2','S50'),
                 ('pos2',int),
                 ('total_copynumber2',int),
                 ('no_chrs_bearing_mutation2',int)]

    mlcn_dtype =[('chr1','S50'),
                 ('pos1',int),
                 ('chr2','S50'),
                 ('pos2',int),
                 ('pos1_bb_CN','S50'),
                 ('pos2_bb_CN','S50'),
                 ('most_likely_ref_copynumber',int),
                 ('most_likely_variant_copynumber',int),
                 ('prop_chrs_bearing_mutation',float),
                 ('support',int),
                 ('depth',int),                
                 ('pv',float),
                 ('ploidy_fixed_pv',float)]

    cmem = np.hstack(clus_members)
    cn_vect = np.empty((0,len(cmem)),dtype=cn_dtype)
    mlcn_vect = np.empty((0,len(cmem)),dtype=mlcn_dtype)
    clus_svs = df.loc[cmem].copy()
    
    #sup,dep,cn_r,cn_v,mu_v,sides,av_cov = cluster.get_sv_vals(df,rlen)
    sup,dep,combos,sides,Nvar = cluster.get_sv_vals(df,rlen)
    phis = clus_cert.average_ccf.values

    for idx,sv in clus_svs.iterrows():
        gtype1,gtype2 = sv['gtype1'].split('|'),sv['gtype2'].split('|')
        gtype1,gtype2 = map(methodcaller('split', ','), gtype1),map(methodcaller('split', ','), gtype2)

        #select the first clonal/major fraction copy-number state as the one to output
        maj_cn1,min_cn1 = map(float,gtype1[0])[:2] if gtype1[0][0]!='' else [0.,0.]
        maj_cn2,min_cn2 = map(float,gtype2[0])[:2] if gtype2[0][0]!='' else [0.,0.]
        
        tot_cn1 = maj_cn1+min_cn2
        tot_cn2 = maj_cn2+min_cn2

        bp1_chr = str(sv['bp1_chr'])
        bp1_pos = int(sv['bp1_pos'])
        bp2_chr = str(sv['bp2_chr'])
        bp2_pos = int(sv['bp2_pos'])
       
        side = sides[idx]
        ref_cn, sc_cn, freq = cluster.get_most_likely_cn(cn_r[idx][side],cn_v[idx][side],\
                                                         mu_v[idx][side],sup[idx],dep[idx],phis[idx],pi) 
        pv = cluster.get_pv(phis[idx],ref_cn,sc_cn,freq,pi)
        pv_pf = cluster.get_pv(phis[idx],2.1,2.1,1/2.1,pi)

        cn_new_row = np.array([(bp1_chr,bp1_pos,tot_cn1,int(sc_cn*freq),bp2_chr,bp2_pos,tot_cn2,int(sc_cn*freq))],dtype=cn_dtype)
        cn_vect = np.append(cn_vect,cn_new_row)
        
        ml_new_row = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos,sv['gtype1'],sv['gtype2'],ref_cn,sc_cn,freq,sup[idx],dep[idx],pv,pv_pf)],dtype=mlcn_dtype)
        mlcn_vect = np.append(mlcn_vect,ml_new_row)
        
    pd.DataFrame(mlcn_vect).to_csv('%s/%s_most_likely_copynumbers.txt'%(clus_out_dir,sample),sep='\t',index=False)
    pd.DataFrame(cn_vect).to_csv('%s/%s_copynumber.txt'%(clus_out_dir,sample),sep='\t',index=False)
    df_probs.to_csv('%s/%s_assignment_probability_table.txt'%(clus_out_dir,sample),sep='\t',index=False)
    clus_cert.to_csv('%s/%s_cluster_certainty.txt'%(clus_out_dir,sample),sep='\t',index=False)

