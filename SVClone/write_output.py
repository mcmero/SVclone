import numpy as np
import pandas as pd
import os
import ipdb
from operator import methodcaller

from . import dtypes
from . import cluster
from . import load_data

def dump_trace(clus_info,center_trace,outf):    
    traces = np.array([center_trace[:,cid] for cid in clus_info.clus_id.values])
    df_traces = pd.DataFrame(np.transpose(traces),columns=clus_info.clus_id)
    df_traces.to_csv(outf,sep='\t',index=False)

def write_out_files(df,clus_info,clus_members,df_probs,clus_cert,clus_out_dir,sample,pi,sup,dep,cn_states,are_snvs=False):
    if are_snvs:
        clus_out_dir = '%s/snvs'%clus_out_dir
        if not os.path.exists(clus_out_dir):
            os.makedirs(clus_out_dir)
    
    clus_info['CCF'] = clus_info.phi.values
    clus_info['phi'] = clus_info.phi.values*pi #convert to proportions
    clus_info = clus_info[['clus_id','size','phi','CCF']]
    rename_cols =  {'clus_id': 'cluster', 'size': 'n_ssms', 'phi': 'proportion'}
    
    clus_info = clus_info.rename(columns = rename_cols)    
    clus_info.to_csv('%s/%s_subclonal_structure.txt'%(clus_out_dir,sample),sep='\t',index=False)

    with open('%s/number_of_clusters.txt'%clus_out_dir,'w') as outf:
        outf.write("sample\tclusters\n")
        outf.write('%s\t%d\n'%(sample,len(clus_info)))

    cmem        = np.hstack(clus_members)
    cn_dtype    = dtypes.sv_cn_dtype    if not are_snvs else dtypes.snv_cn_dtype
    mlcn_dtype  = dtypes.sv_mlcn_dtype  if not are_snvs else dtypes.snv_mlcn_dtype
    mult_dtype  = dtypes.sv_mult_dtype  if not are_snvs else dtypes.snv_mult_dtype

    cn_vect     = np.empty((0,len(cmem)),dtype=cn_dtype)
    mlcn_vect   = np.empty((0,len(cmem)),dtype=mlcn_dtype)
    mult_vect   = np.empty((0,len(cmem)),dtype=mult_dtype)
    clus_vars    = df.loc[cmem].copy()
    
    phis = clus_cert.average_ccf.values
    cns,pvs = cluster.get_most_likely_cn_states(cn_states,sup,dep,phis,pi)

    for idx,var in clus_vars.iterrows():
        gtype1,gtype2 = None,None
        if not are_snvs:
            gtype1,gtype2 = var['gtype1'].split('|'),var['gtype2'].split('|')
            gtype1,gtype2 = map(methodcaller('split', ','), gtype1),map(methodcaller('split', ','), gtype2)

            #select the first clonal/major fraction copy-number state as the one to output
            maj_cn1,min_cn1 = map(float,gtype1[0])[:2] if gtype1[0][0]!='' else [0.,0.]
            maj_cn2,min_cn2 = map(float,gtype2[0])[:2] if gtype2[0][0]!='' else [0.,0.]

            maj_cn1,min_cn1 = map(float,gtype1[0])[:2] if gtype1[0][0]!='' else [0.,0.]
            tot_cn1 = maj_cn1 + min_cn1
            tot_cn2 = maj_cn2+min_cn2 if not are_snvs else 0

            bp1_chr = str(var['bp1_chr'])
            bp1_pos = int(var['bp1_pos'])
            bp2_chr = str(var['bp2_chr'])
            bp2_pos = int(var['bp2_pos'])

            ref_cn, sc_cn, freq = cns[idx]
            pv = pvs[idx]
            chrs_bearing_mut = int(sc_cn*freq) if (not np.isnan(sc_cn) or not np.isnan(freq)) else 0

            ref_cn = 0 if np.isnan(sc_cn) else ref_cn
            sc_cn = 0 if np.isnan(sc_cn) else sc_cn
            freq = 0 if np.isnan(freq) else freq

            cn_new_row = np.array([(bp1_chr,bp1_pos,tot_cn1,chrs_bearing_mut,
                                    bp2_chr,bp2_pos,tot_cn2,chrs_bearing_mut)],dtype=cn_dtype)
            ml_new_row = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos,var['gtype1'],var['gtype2'],
                                    ref_cn,sc_cn,freq,sup[idx],dep[idx],pv)],dtype=mlcn_dtype)
            
            var_states = cn_states[idx]
            tot_opts = ','.join(map(str,[int(x[1]) for x in var_states]))
            var_opts = ','.join(map(str,[int(round(x[1]*x[2])) for x in var_states]))
            probs =  cluster.get_probs(var_states,sup[idx],dep[idx],phis[idx],pi)
            
            mult_new_row = np.array([(bp1_chr,bp1_pos,bp2_chr,bp2_pos,sc_cn,
                                     int(round(freq*sc_cn)),tot_opts,var_opts,probs)],dtype=mult_dtype)
            
            cn_vect   = np.append(cn_vect,cn_new_row)
            mlcn_vect = np.append(mlcn_vect,ml_new_row)            
            mult_vect = np.append(mult_vect,mult_new_row)
        else:
            gtype1 = map(lambda x: x.split(','), var['gtype'].split('|'))

            maj_cn1,min_cn1 = map(float,gtype1[0])[:2] if gtype1[0][0]!='' else [0.,0.]
            tot_cn1 = maj_cn1 + min_cn1

            chrom = str(var['chrom'])
            pos = int(var['pos'])
       
            ref_cn, sc_cn, freq = cns[idx]
            pv = pvs[idx]

            cn_new_row = np.array([(chrom,pos,tot_cn1,int(sc_cn*freq))],dtype=cn_dtype)
            ml_new_row = np.array([(chrom,pos,var['gtype'],ref_cn,sc_cn,freq)],dtype=mlcn_dtype)

            var_states = cn_states[idx]
            tot_opts = ','.join(map(str,[int(x[1]) for x in var_states]))
            var_opts = ','.join(map(str,[int(round(x[1]*x[2])) for x in var_states]))
            probs =  cluster.get_probs(var_states,sup[idx],dep[idx],phis[idx],pi)

            mult_new_row = np.array([(chrom,pos,sc_cn,int(round(freq*sc_cn)),tot_opts,var_opts,probs)],dtype=mult_dtype)
            
            cn_vect = np.append(cn_vect,cn_new_row)    
            mlcn_vect = np.append(mlcn_vect,ml_new_row)
            mult_vect = np.append(mult_vect,mult_new_row)

    #adjust cluster freqs to cell prevalence
    clus_cert.average_ccf = clus_cert.average_ccf.values*pi
    clus_cert['90_perc_CI_lo'] = clus_cert['90_perc_CI_lo'].values*pi
    clus_cert['90_perc_CI_hi'] = clus_cert['90_perc_CI_hi'].values*pi

    #rename cols
    rename_cols =  {'bp1_chr': 'chr1', 'bp1_pos': 'pos1', 
                    'bp2_chr': 'chr2', 'bp2_pos': 'pos2', 
                    'average_ccf': 'average_proportion'}
    clus_cert = clus_cert.rename(columns = rename_cols)
    df_probs = df_probs.rename(columns = rename_cols)

    pd.DataFrame(mlcn_vect).to_csv('%s/%s_most_likely_copynumbers.txt'%(clus_out_dir,sample),sep='\t',index=False)
    pd.DataFrame(mult_vect).to_csv('%s/%s_multiplicity.txt'%(clus_out_dir,sample),sep='\t',index=False)
    pd.DataFrame(cn_vect).to_csv('%s/%s_copynumber.txt'%(clus_out_dir,sample),sep='\t',index=False)
    df_probs.to_csv('%s/%s_assignment_probability_table.txt'%(clus_out_dir,sample),sep='\t',index=False)
    clus_cert.to_csv('%s/%s_cluster_certainty.txt'%(clus_out_dir,sample),sep='\t',index=False)

