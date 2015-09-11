import numpy as np
import pandas as pd
import os
import ipdb
import math
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
    
    clus_info['phi'] = clus_info.phi.values*pi
    clus_info = clus_info[['clus_id','size','phi']]
    rename_cols =  {'clus_id': 'cluster', 'size': 'n_ssms', 'phi': 'proportion'}

    clus_info = clus_info.rename(columns = rename_cols)    
    clus_info.to_csv('%s/%s_subclonal_structure.txt'%(clus_out_dir,sample),sep='\t',index=False)

    with open('%s/number_of_clusters.txt'%clus_out_dir,'w') as outf:
        outf.write("sample\tclusters\n")
        outf.write('%s\t%d\n'%(sample,len(clus_info)))
    
    cn_dtype =  [('chr','S50'),
                 ('pos',int),
                 ('total_copynumber',int),
                 ('no_chrs_bearing_mutation',int)]

    mlcn_dtype =[('chr','S50'),
                 ('pos',int),
                 ('bb_CN','S50'),
                 ('most_likely_ref_copynumber',int),
                 ('most_likely_variant_copynumber',int),
                 ('prop_chrs_bearing_mutation',float)]
    
    mult_dtype =[('chr','S50'),
                 ('pos',int),
                 ('tumour_copynumber',int),
                 ('multiplicity',int),
                 ('tumour_copynumber_options','S50'),
                 ('multiplicity_options','S50'),
                 ('probabilities','S50')]

    cmem = np.hstack(clus_members)
    cn_vect = np.empty((0,len(cmem)),dtype=cn_dtype)
    mlcn_vect = np.empty((0,len(cmem)),dtype=mlcn_dtype)
    mult_vect   = np.empty((0,len(cmem)),dtype=mult_dtype)
    #clus_snvs = df.loc[cmem].copy()
    df = df.fillna('')
    df = df[df.chrom.values!='']

    #sup,n,cn_r,cn_v,mu_v = cluster.get_snv_vals(df)
    sup,dep,cn_states,Nvar = cluster.get_snv_vals(df)
    phis = clus_cert.average_ccf.values
    cns, pvs = cluster.get_most_likely_cn_states(cn_states,sup,dep,phis,pi)

    for idx,snv in df.iterrows():
        gtype = snv['gtype'].split('|')
        gtype = map(methodcaller('split', ','), gtype)

        maj_cn,min_cn = map(float,gtype[0])[:2] if gtype[0][0]!='' else [0.,0.]
        tot_cn = maj_cn+min_cn

        chrom = str(snv['chrom'])
        pos = int(snv['pos'])
        #ref_cn, sc_cn, freq = cluster.get_most_likely_cn(cn_r[idx][0],cn_v[idx][0],\
        #                                                 mu_v[idx][0],sup[idx],dep[idx],phis[idx],pi) 
        ref_cn, sc_cn, freq = cns[idx]
        pv = pvs[idx]
        
        cn_new_row = np.array([(chrom,pos,tot_cn,int(sc_cn*freq))],dtype=cn_dtype)
        cn_vect = np.append(cn_vect,cn_new_row)
        
        ml_new_row = np.array([(chrom,pos,snv['gtype'],ref_cn,sc_cn,freq)],dtype=mlcn_dtype)
        mlcn_vect = np.append(mlcn_vect,ml_new_row)
        
        var_states = cn_states[idx]
        tot_opts = ','.join(map(str,[int(x[1]) for x in var_states]))
        var_opts = ','.join(map(str,[int(round(x[1]*x[2])) for x in var_states]))
        probs = [(1/float(len(var_states)))]*len(var_states)
        probs = ','.join(map(lambda x: str(round(x,4)),probs))

        mult_new_row = np.array([(chrom,pos,sc_cn,int(round(freq*sc_cn)),tot_opts,var_opts,probs)],dtype=mult_dtype)
        mult_vect = np.append(mult_vect,mult_new_row)
        
    clus_cert.average_ccf = clus_cert.average_ccf.values*pi
    clus_cert['90_perc_CI_lo'] = clus_cert['90_perc_CI_lo'].values*pi
    clus_cert['90_perc_CI_hi'] = clus_cert['90_perc_CI_hi'].values*pi

    #rename cols
    rename_cols =  {'bp1_chr': 'chr1', 'bp1_pos': 'pos1', 'bp2_chr': 'chr2', 
                    'bp2_pos': 'pos2', 'average_ccf': 'average_proportion'}
    clus_cert = clus_cert.rename(columns = rename_cols)
    df_probs = df_probs.rename(columns = rename_cols)

    pd.DataFrame(mlcn_vect).to_csv('%s/%s_most_likely_copynumbers.txt'%(clus_out_dir,sample),sep='\t',index=False)
    pd.DataFrame(cn_vect).to_csv('%s/%s_copynumber.txt'%(clus_out_dir,sample),sep='\t',index=False)
    df_probs.to_csv('%s/%s_assignment_probability_table.txt'%(clus_out_dir,sample),sep='\t',index=False)
    clus_cert.to_csv('%s/%s_cluster_certainty.txt'%(clus_out_dir,sample),sep='\t',index=False)
    pd.DataFrame(mult_vect).to_csv('%s/%s_multiplicity.txt'%(clus_out_dir,sample),sep='\t',index=False)

def write_out_files(df,clus_info,clus_members,df_probs,clus_cert,clus_out_dir,sample,pi,ploidy,rlen):
    
    clus_info['phi'] = clus_info.phi.values*pi
    clus_info = clus_info[['clus_id','size','phi']]
    rename_cols =  {'clus_id': 'cluster', 'size': 'n_ssms', 'phi': 'proportion'}

    clus_info = clus_info.rename(columns = rename_cols)    
    clus_info.to_csv('%s/%s_subclonal_structure.txt'%(clus_out_dir,sample),sep='\t',index=False)

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
                 ('pv',float)]
                 #('ploidy_fixed_pv',float)]
   
    mult_dtype =[('chr1','S50'),
                 ('pos1',int),
                 ('chr2','S50'),
                 ('pos2',int),
                 ('tumour_copynumber',int),
                 ('multiplicity',int),
                 ('tumour_copynumber_options','S50'),
                 ('multiplicity_options','S50'),
                 ('probabilities','S50')]

    cmem        = np.hstack(clus_members)
    cn_vect     = np.empty((0,len(cmem)),dtype=cn_dtype)
    mlcn_vect   = np.empty((0,len(cmem)),dtype=mlcn_dtype)
    mult_vect   = np.empty((0,len(cmem)),dtype=mult_dtype)
    clus_svs    = df.loc[cmem].copy()
    
    #sup,dep,cn_r,cn_v,mu_v,sides,av_cov = cluster.get_sv_vals(df,rlen)
    sup,dep,cn_states,Nvar = cluster.get_sv_vals(df,rlen,pi,ploidy)
    #TODO: should all be cellular prevalence - do this with SNVs also 
    phis = clus_cert.average_ccf.values
    cns,pvs = cluster.get_most_likely_cn_states(cn_states,sup,dep,phis,pi)

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
       
        #side = sides[idx]
        #ref_cn, sc_cn, freq = cluster.get_most_likely_cn(cn_r[idx][side],cn_v[idx][side],\
        #                                                 mu_v[idx][side],sup[idx],dep[idx],phis[idx],pi) 
        #pv = cluster.get_pv(phis[idx],ref_cn,sc_cn,freq,pi)
        #pv_pf = cluster.get_pv(phis[idx],2.1,2.1,1/2.1,pi)
        ref_cn, sc_cn, freq = cns[idx]
        pv = pvs[idx]

        cn_new_row = np.array([(bp1_chr,bp1_pos,tot_cn1,int(sc_cn*freq),
                                bp2_chr,bp2_pos,tot_cn2,int(sc_cn*freq))],dtype=cn_dtype)
        cn_vect = np.append(cn_vect,cn_new_row)
        
        ml_new_row = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos,sv['gtype1'],sv['gtype2'],
                                ref_cn,sc_cn,freq,sup[idx],dep[idx],pv)],dtype=mlcn_dtype)
        mlcn_vect = np.append(mlcn_vect,ml_new_row)

        var_states = cn_states[idx]
        tot_opts = ','.join(map(str,[int(x[1]) for x in var_states]))
        var_opts = ','.join(map(str,[int(round(x[1]*x[2])) for x in var_states]))
        #TODO: convert to true probabilities
        #probs = cluster.calc_lik(var_states,sup[idx],dep[idx],phis[idx],pi)[1]
        #probs = map(math.exp,probs)
        #probs = np.array(probs)/sum(probs)
        probs = [(1/float(len(var_states)))]*len(var_states)
        probs = ','.join(map(lambda x: str(round(x,4)),probs))
        mult_new_row = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos,sc_cn,
                                  int(round(freq*sc_cn)),tot_opts,var_opts,probs)],dtype=mult_dtype)
        mult_vect = np.append(mult_vect,mult_new_row)
    
    #adjust cluster freqs to cell prevalence
    clus_cert.average_ccf = clus_cert.average_ccf.values*pi
    clus_cert['90_perc_CI_lo'] = clus_cert['90_perc_CI_lo'].values*pi
    clus_cert['90_perc_CI_hi'] = clus_cert['90_perc_CI_hi'].values*pi

    #rename cols
    rename_cols =  {'bp1_chr': 'chr1', 'bp1_pos': 'pos1', 'bp2_chr': 'chr2', 'bp2_pos': 'pos2', 'average_ccf': 'average_proportion'}
    clus_cert = clus_cert.rename(columns = rename_cols)
    df_probs = df_probs.rename(columns = rename_cols)

    pd.DataFrame(mlcn_vect).to_csv('%s/%s_most_likely_copynumbers.txt'%(clus_out_dir,sample),sep='\t',index=False)
    pd.DataFrame(mult_vect).to_csv('%s/%s_multiplicity.txt'%(clus_out_dir,sample),sep='\t',index=False)
    pd.DataFrame(cn_vect).to_csv('%s/%s_copynumber.txt'%(clus_out_dir,sample),sep='\t',index=False)
    df_probs.to_csv('%s/%s_assignment_probability_table.txt'%(clus_out_dir,sample),sep='\t',index=False)
    clus_cert.to_csv('%s/%s_cluster_certainty.txt'%(clus_out_dir,sample),sep='\t',index=False)

