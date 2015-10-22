'''
Run clustering and tree building on sample inputs
'''
from __future__ import print_function

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
from . import load_data
from . import parameters as params

def get_outlier_ranges(vals):
    q1 = np.percentile(vals,25)
    q3 = np.percentile(vals,75)
    iqr = max(q3 - q1, 1)
    lower = q1 - 1.5*iqr
    upper = q3 + 1.5*iqr
    return [lower,upper]

def run_simple_filter(df,rlen,insert,minsplit,minspan,sizefilter,min_dep,valid_chrs):
    '''
    filter based on presence of supporting reads and > fragment length
    ''' 
    span = np.array(df.spanning)
    split = np.array(df.bp1_split+df.bp2_split)
    #df_flt = df[(span+split)>=1]
    df_flt = df[np.logical_and(span>=minspan,split>=minsplit)]

    dep1 = df_flt.support + df_flt.norm1
    dep2 = df_flt.support + df_flt.norm2
    df_flt = df_flt[np.logical_and(dep1>=min_dep,dep2>=min_dep)]

    print('Filtered out %d SVs based on minimum depth and spanning/split read limits' % (len(df) - len(df_flt)))

    itx = df_flt['bp1_chr']!=df_flt['bp2_chr']
    frag_len = 2*rlen+insert
    sizefilter = sizefilter if sizefilter>=0 else frag_len
    df_flt2 = df_flt[ itx | (abs(df_flt['bp1_pos']-df_flt['bp2_pos'])>sizefilter) ]

    print('Filtered out %d SVs based on size limits' % (len(df_flt) - len(df_flt2)))
    
    if valid_chrs:
        chr1 = np.array([chrom in params.valid_chroms for chrom in df_flt2.bp1_chr])
        chr2 = np.array([chrom in params.valid_chroms for chrom in df_flt2.bp2_chr])
        df_flt3 = df_flt2[np.logical_and(chr1,chr2)]
        
        print('Filtered out %d SVs that had non-standard chromosomes' % (len(df_flt2) - len(df_flt3)))
        df_flt2 = pd.DataFrame(df_flt3,copy=True)

    return df_flt2

def is_clonal_neutral(gtype):
    if gtype=='': return False
    gtype = map(methodcaller('split',','),gtype.split('|'))
    if len(gtype)==1:
        gt = map(float,gtype[0])
        return (gt[0]==1 and gt[1]==1 and gt[2]==1)
    return False

def is_copynumber_zero(gtype):
    '''
    returns false if copynumber genotype
    of the major and minor alleles
    are 0 for the clone or any subclone
    '''
    if gtype=='': return False
    gtype = map(methodcaller('split',','),gtype.split('|'))
    if len(gtype)==1:
        gt = map(float,gtype[0])
        if gt[2]==1:            
            return (gt[0]==0 and gt[1]==0)
        else:
            gt2 = map(float,gtype[1])
            return ((gt[0]==0 and gt[1]==0) or (gt2[0]==0 and gt2[1]==0))
    return False

def get_weighted_cns(gtypes):
    gtypes_split = [map(methodcaller("split",","),x) for x in map(methodcaller("split","|"),gtypes)]    
    cn_vals = []
    for gtype in gtypes_split:
        cn_val = sum([int(eval(g[0]))+int(eval(g[1]))*float(eval(g[2])) if g!=[''] else 2 for g in gtype])
        cn_vals.append(cn_val)
    return np.array(cn_vals)/2

def normalise_wins_by_cn(df_flt):
    bp1_win = df_flt.bp1_win_norm.map(float).values
    bp2_win = df_flt.bp2_win_norm.map(float).values
    
    bp1_wcn = get_weighted_cns(df_flt.gtype1.values)
    bp2_wcn = get_weighted_cns(df_flt.gtype2.values)

    bp1_nonzero = np.logical_not(bp1_wcn==0)
    bp2_nonzero = np.logical_not(bp2_wcn==0)
    bp1_win[bp1_nonzero] = (bp1_win[bp1_nonzero]/bp1_wcn[bp1_nonzero])
    bp2_win[bp2_nonzero] = (bp2_win[bp2_nonzero]/bp2_wcn[bp2_nonzero])
    
    return bp1_win,bp2_win

def filter_outlying_norm_wins(df_flt):
    bp1_win, bp2_win = normalise_wins_by_cn(df_flt)
    bp1_ranges = get_outlier_ranges(bp1_win)
    bp2_ranges = get_outlier_ranges(bp2_win)

    bp1_flt = np.logical_and(bp1_win>bp1_ranges[0],bp1_win<bp1_ranges[1])
    bp2_flt = np.logical_and(bp2_win>bp2_ranges[0],bp2_win<bp2_ranges[1])
    df_flt = df_flt[np.logical_and(bp1_flt,bp2_flt)]
    
    return df_flt

def run_cnv_filter(df_flt,cnv,neutral,filter_outliers,are_snvs=False):
    '''
    filter based on either CNV neutral, or presence of CNV vals
    '''
    n_df = len(df_flt)
    if len(cnv)>0 and neutral:
        # filter out copy-aberrant SVs and outying norm read counts (>1-percentile)
        # major and minor copy-numbers must be 1
        is_neutral = []
        if are_snvs:
            df_flt = df_flt.fillna('')
            df_flt = df_flt[df_flt.gtype.values!='']
            is_neutral = map(is_clonal_neutral,df_flt.gtype.values)
            df_flt2 = df_flt[is_neutral]
            print('Filtered out %d SNVs that were not copy-number neutral' % (n_df - len(df_flt)))

            if filter_outliers:
                n_df = len(df_flt)
                depths = df_flt['ref'].values+df_flt['var'].values
                dep_ranges = get_outlier_ranges(depths)
                df_flt = df_flt[np.logical_and(depths>dep_ranges[0],depths<dep_ranges[1])]
                print('Filtered out %d SNVs which had outlying depths' % (n_df - len(df_flt)))
        else: 
            gt1_is_neutral = map(is_clonal_neutral,df_flt.gtype1.values)
            gt2_is_neutral = map(is_clonal_neutral,df_flt.gtype2.values)
            is_neutral = np.logical_and(gt1_is_neutral,gt2_is_neutral)
            df_flt = df_flt[is_neutral] 
            print('Filtered out %d SVs which were not copy-number neutral' % (n_df - len(df_flt)))

            if filter_outliers:
                n_df = len(df_flt)
                df_flt = filter_outlying_norm_wins(df_flt)
                print('Filtered out %d SVs which had outlying depths' % (n_df - len(df_flt)))

    elif len(cnv)>0:
        if are_snvs:
            df_flt = df_flt.fillna('')
            df_flt = df_flt[df_flt.gtype.values!='']
            print('Filtered out %d SNVs which did not have copy-numbers' % (n_df - len(df_flt)))

            if filter_outliers:
                # weight ranges by copy-numbers
                depths = df_flt['ref'].values + df_flt['var'].values
                cns = get_weighted_cns(df_flt.gtype.values)
                cn_nonzero = np.logical_not(cns==0)
                depths[cn_nonzero] = (depths[cn_nonzero]/cns[cn_nonzero])
            
                n_df = len(df_flt)
                dep_ranges = get_outlier_ranges(depths)
                df_flt = df_flt[np.logical_and(depths>dep_ranges[0],depths<dep_ranges[1])]
                print('Filtered out %d SNVs which had outlying depths' % (n_df - len(df_flt)))
        else:
            df_flt = df_flt[np.logical_or(df_flt.gtype1.values!='',df_flt.gtype2.values!='')]
            
            gt1_is_zero = np.array(map(is_copynumber_zero,df_flt.gtype1.values))
            gt2_is_zero = np.array(map(is_copynumber_zero,df_flt.gtype2.values))
            df_flt = df_flt[np.logical_and(gt1_is_zero==False,gt2_is_zero==False)]
            print('Filtered out %d SVs without copy-numbers or with 0,0 copy-numbers' % (n_df - len(df_flt)))
            n_df = len(df_flt)            

            if filter_outliers:
                df_flt = filter_outlying_norm_wins(df_flt)
                print('Filtered out %d SVs which had outlying depths' % (n_df - len(df_flt)))

    return df_flt

def match_copy_numbers(var_df, cnv_df, bp_fields=['bp1_chr','bp1_pos'], gtype_field='gtype1'):
    
    var_df[gtype_field] = ''
    chrom_field, pos_field = bp_fields
    var_df[chrom_field] = map(str,var_df[chrom_field].values)
    var_df[pos_field] = map(int,var_df[pos_field].values)
    bp_chroms = np.unique(var_df[chrom_field].values)
    bp_chroms = sorted(bp_chroms, key=lambda item: (int(item.partition(' ')[0]) \
                        if item[0].isdigit() else float('inf'), item))
    
    for bchr in bp_chroms:
        #print('Matching copy-numbers for chrom %s'%bchr)
        gtypes = []
        current_chr = var_df[chrom_field].values==bchr
        var_tmp = var_df[current_chr]
        cnv_tmp = cnv_df[cnv_df['chr']==bchr]
        
        if len(cnv_tmp)==0:
            continue

        for pos in var_tmp[pos_field].values:
            #print("Checking overlaps for pos %d\n" %pos)
            cnv_start_list = cnv_tmp.startpos.values
            cnv_end_list   = cnv_tmp.endpos.values
            overlaps = np.logical_and(pos >= cnv_start_list, pos <= cnv_end_list)

            match = cnv_tmp[overlaps]            
            if len(match)==0:
                # TODO: matching of SV/CNV boundaries
                # right now just takes the position that's 
                # the closest CNV value less than the SV 
                # position or the first more than if there 
                # are no CNV segments that are less than
                if np.any(pos >= cnv_end_list):
                    greater_than = cnv_tmp[pos >= cnv_end_list].index
                    closest = greater_than[len(greater_than)-1]                    
                    gtypes.append(cnv_tmp.loc[closest].gtype)
                elif np.any(pos <= cnv_start_list):
                    less_than = cnv_tmp[pos <= cnv_start_list].index
                    closest = less_than[0]
                    gtypes.append(cnv_tmp.loc[closest].gtype)
                else:
                    gtypes.append('')
            else:
                gtype = match.loc[match.index[0]].gtype
                gtypes.append(gtype)
        
        var_indexes = var_df[current_chr].index.values
        var_df.loc[var_indexes,gtype_field] = gtypes
    return var_df

def is_same_sv_germline(sv1,sv2):
    sv1_chr1, sv1_bp1, sv1_chr2, sv1_bp2 = sv1
    sv2_chr1, sv2_bp1, sv2_chr2, sv2_bp2 = sv2

    if sv1_chr1==sv2_chr1 and sv1_chr2==sv2_chr2:
        if abs(sv1_bp1-sv2_bp1)<params.gl_th and abs(sv1_bp2-sv2_bp2)<params.gl_th:
            return True
    if sv1_chr2==sv2_chr1 and sv1_chr1==sv2_chr2:
        if abs(sv1_bp2-sv2_bp1)<params.gl_th and abs(sv1_bp1-sv2_bp2)<params.gl_th:
            return True
    return False

def filter_germline(gml_file,sv_df,rlen,insert):
    print("Filtering out germline SVs...")
    df_gml = pd.DataFrame(pd.read_csv(gml_file,delimiter='\t',dtype=None,low_memory=False))
    #df_gml = run_simple_filter(df_gml,rlen,insert)
    germline = []
    
    for idx_sv,sv in sv_df.iterrows():
        sv_loc = [str(sv.bp1_chr),int(sv.bp1_pos),str(sv.bp2_chr),int(sv.bp2_pos)]
        same_chrs = np.logical_and(df_gml.bp1_chr.values==sv_loc[0],df_gml.bp2_chr.values==sv_loc[2])
        df_gml_tmp = df_gml[same_chrs]
        for idx_gml,sv_gml in df_gml_tmp.iterrows():
            sv_loc_gml = [str(sv_gml.bp1_chr),int(sv_gml.bp1_pos),str(sv_gml.bp2_chr),int(sv_gml.bp2_pos)]
            if sv_gml.support>0 and is_same_sv_germline(sv_loc,sv_loc_gml):
                germline.append(idx_sv)
                break

    print("Filtered out %d SVs that were found in the germline!" % len(germline))
    return sv_df.drop(germline,axis=0)


def adjust_sv_read_counts(sv_df,pi,pl,min_dep):
    n = zip(np.array(sv_df.norm1.values),np.array(sv_df.norm2.values))
    s = np.array(sv_df.bp1_split.values+sv_df.bp2_split.values)
    d = np.array(sv_df.spanning.values)
    sup = np.array(d+s,dtype=float)
    Nvar = len(sv_df)
    
    # if low normal depth on one side only, pick other side 
    sides = np.zeros(Nvar,dtype=int)
    low_norm_dep = np.logical_xor(sv_df.norm1.values<min_dep,sv_df.norm2.values<min_dep)
    low_norm_vals = zip(sv_df.norm1.values[low_norm_dep],sv_df.norm2.values[low_norm_dep])
    sides[low_norm_dep] = [ 0 if n1 < min_dep else 1 for n1, n2 in low_norm_vals ]
    
    # if one side doesn't have CNV data, pick the other side
    sides[sv_df.gtype1.values==''] = 1 
    sides[sv_df.gtype2.values==''] = 0

    # prefer sides with subclonal genotype data    
    gt1_sc = np.array(map(len,map(methodcaller("split","|"),sv_df.gtype1.values)))>1
    gt2_sc = np.array(map(len,map(methodcaller("split","|"),sv_df.gtype1.values)))>1    
    one_sc = np.logical_xor(gt1_sc,gt2_sc)

    combos = sv_df.apply(cluster.get_sv_allele_combos,axis=1)
    exclusive_subclones = zip(sv_df.gtype1.values[one_sc],sv_df.gtype2.values[one_sc]) 
    sides[one_sc] = [0 if gt1!='' else 1 for gt1,gt2 in exclusive_subclones]
    has_both_gts = np.logical_and(sv_df.gtype1.values!='',sv_df.gtype2.values!='')
    cn_states = [cn[side] for cn,side in zip(combos,sides)]
    
    # both sides have genotypes, both either subclonal or clonal
    # in this case, just take the simple normal mean of the two sides
    norm = np.array([float(ni[si]) for ni,si in zip(n,sides)])
    both_gts_same_type = np.logical_and(has_both_gts,one_sc==False)
    norm[both_gts_same_type] = map(np.mean,np.array(n)[both_gts_same_type])

    try:

        # read value adjustments for specific types of events
        # currently the adjustments are quite simple
        sv_classes = sv_df.classification.values
        #invs = [ idx in params.inversion_class for idx,sv_class in enumerate(sv_classes) ]
        dups = np.array([ sv_class in params.dna_gain_class for idx,sv_class in enumerate(sv_classes) ])
        dels = np.array([ sv_class in params.deletion_class for idx,sv_class in enumerate(sv_classes) ])
        
        # normal read counts for duplications are adjusted by purity and ploidy
        if sum(dups)>0:
            # estimate adjustment from normal counts for SV events where there is no gain of DNA
            # if these events don't exist, adjust by normal component + tumour/2 (half tumour normal
            # reads will come from duplication, on one chromosome, so divide by ploidy
            adjust_factor = np.mean(norm[dels])/np.mean(norm[dups]) if sum(dels)>0 else (1-float(pi)) + (float(pi)/2/pl)
            norm[dups]    = norm[dups] * adjust_factor
        
    except AttributeError:
        print('Warning, no valid classifications found. SV read counts cannot be adjusted')

    sv_df['adjusted_norm']      = norm
    sv_df['adjusted_support']   = sup
    sv_df['adjusted_depth']     = norm+sup
    sv_df['preferred_side']     = sides
    sv_df['raw_mean_vaf']       = (s+d)/(s+d+map(float,sv_df.norm_mean.values))
    sv_df['adjusted_vaf']       = map(float,sup)/(norm+sup)

    return sv_df

def run(args):    
    # pr = cProfile.Profile()
    # pr.enable()    

    sample      = args.sample
    svs         = args.procd_svs
    gml         = args.germline
    cnvs        = args.cnvs
    out         = args.outdir
    params_file = args.params_file 
    snvs        = args.snvs
    snv_format  = args.snv_format
    neutral     = args.neutral
    pi          = args.pi
    ploidy      = args.ploidy
    minsplit    = int(args.minsplit)
    minspan     = int(args.minspan)
    sizefilter  = int(args.sizefilter)
    filter_otl  = args.filter_outliers
    min_dep     = args.min_dep
    valid_chrs  = args.valid_chrs
   
    def proc_arg(arg,n_args=1,of_type=str):
        arg = str.split(arg,',')
        arg = arg * n_args if len(arg)==1 else arg
        if of_type==int or of_type==float:
            return map(eval,arg)
        else:
            return map(of_type,arg)

    sample  = proc_arg(sample)
    n       = len(sample)
    svs     = proc_arg(svs)
    cnvs    = proc_arg(cnvs)
    pi      = proc_arg(pi,n,float)
    ploidy  = proc_arg(ploidy,n,float)    

    rlen    = 100
    insert  = 100
    
    for p in pi:
        if p<0 or p>1:
            raise ValueError("Tumour purity value not between 0 and 1!")
            
    if len(svs)!=n or len(cnvs)!=n:
        raise ValueError("Number of samples does not match number of input files")

    if not os.path.exists(out):
        os.makedirs(out)

    if len(sample)>1:
        print("Multiple sample processing not yet implemented")
        #TODO: processing of multiple samples
    else:
        sample,sv,cnv,pi,ploidy = sample[0],svs[0],cnvs[0],pi[0],ploidy[0]
        sv_df  = pd.DataFrame()
        cnv_df = pd.DataFrame()
        snv_df = pd.DataFrame()
       
        if params_file=='':
            params_file = '%s/%s_params.txt' % (out,sample) 
        if os.path.exists(params_file):
            read_params = pd.read_csv(params_file,delimiter='\t',dtype=None,header=0)
            rlen        = int(read_params.read_len.values[0])
            insert      = float(read_params.insert_mean.values[0])
        else:
            print('%s/%s_params.txt file not found! Assuming read length = 100, mean insert length = 100' % (out,sample)) 

        if snvs!="":
            if snv_format == 'sanger':
                snv_df = load_data.load_snvs_sanger(snvs)
            elif snv_format == 'mutect':
                snv_df = load_data.load_snvs_mutect(snvs,sample)
            elif snv_format == 'mutect_callstats':
                snv_df = load_data.load_snvs_mutect_callstats(snvs)
            dep = snv_df['ref'] + snv_df['var']            
            snv_df = snv_df[dep>=min_dep]
    
        if sv!="":
            sv_df = load_data.load_svs(sv)
            sv_df = run_simple_filter(sv_df,rlen,insert,minsplit,minspan,sizefilter,min_dep,valid_chrs)
        
        if gml!="":
            sv_df = filter_germline(gml,sv_df,rlen,insert)
       
        if cnv!="":
            cnv_df = load_data.load_cnvs(cnv)

            if len(sv_df)>0:
                print('Matching copy-numbers for SVs...')
                sv_df = match_copy_numbers(sv_df,cnv_df) 
                sv_df = match_copy_numbers(sv_df,cnv_df,['bp2_chr','bp2_pos'],'gtype2') 
                sv_df = run_cnv_filter(sv_df,cnv,neutral,filter_otl)
                print('Keeping %d SVs' % len(sv_df))
        
            if len(snv_df)>0:
                print('Matching copy-numbers for SNVs...')
                snv_df = match_copy_numbers(snv_df,cnv_df,['chrom','pos'],'gtype')
                snv_df = run_cnv_filter(snv_df,cnv,neutral,filter_otl,are_snvs=True)
                print('Keeping %d SNVs' % len(snv_df))
        else:

            if len(sv_df)>0:
                print('No Battenberg input defined, assuming all loci are copy-number neutral')
                sv_df['gtype1'] = '1,1,1.000000'
                sv_df['gtype2'] = '1,1,1.000000'
                sv_df = run_cnv_filter(sv_df,cnv,neutral,filter_otl)
                print('Retained %d SVs' % len(sv_df))

            if len(snv_df)>0:
                snv_df['gtype'] = '1,1,1.000000'
                snv_df = run_cnv_filter(snv_df,cnv,neutral,filter_otl,are_snvs=True)
                print('Retained %d SNVs' % len(snv_df))
            
        if len(sv_df)>0: 
            sv_df.index = range(len(sv_df)) #reindex
            sv_df = adjust_sv_read_counts(sv_df,pi,ploidy,min_dep)
            sv_df.to_csv('%s/%s_filtered_svs.tsv'%(out,sample),sep='\t',index=False,na_rep='')
        
        if len(snv_df)>0: 
            snv_df.index = range(len(snv_df)) #reindex
            snv_df.to_csv('%s/%s_filtered_snvs.tsv'%(out,sample),sep='\t',index=False,na_rep='')
        
        with open('%s/purity_ploidy.txt'%out,'w') as outf:
            outf.write("sample\tpurity\tploidy\n")
            outf.write('%s\t%f\t%f\n'%(sample,pi,ploidy))
        
        with open('%s/read_params.txt'%out,'w') as outf:
            outf.write("sample\tread_len\tinsert_mean\n")
            outf.write('%s\t%f\t%f\n'%(sample,rlen,insert))
        
        # pr.disable()
        # s = StringIO.StringIO()
        # sortby = 'cumulative'
        # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        # ps.print_stats()
        # print s.getvalue()

