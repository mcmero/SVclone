'''
Run clustering and tree building on sample inputs
'''
from __future__ import print_function

import os
import configparser
import numpy as np
import re
import pandas as pd
import vcf

from operator import methodcaller
from SVclone.SVprocess import svp_load_data as svp_load
from . import run_clus
from . import load_data

pd.options.mode.chained_assignment = None

def get_outlier_ranges(vals):
    q1 = np.percentile(vals,25)
    q3 = np.percentile(vals,75)
    iqr = max(q3 - q1, 1)
    lower = q1 - 1.5*iqr
    upper = q3 + 1.5*iqr
    return [lower,upper]

def run_simple_filter(df,rlen,insert,minsplit,minspan,sizefilter,min_dep,filter_chrs,valid_chrs,blist):
    '''
    filter based on presence of supporting reads and > fragment length
    '''
    df_flt = df.copy()

    n = len(df_flt)
    itx = df_flt['chr1'] != df_flt['chr2']
    frag_len = 2 * rlen + insert
    sizefilter = sizefilter if sizefilter>=0 else frag_len
    df_flt = df_flt[ itx | (abs(df_flt['pos1'] - df_flt['pos2']) > sizefilter) ]
    print('Filtered out %d SVs based on size limits' % (n - len(df_flt)))

    n = len(df_flt)
    dep1 = df_flt.support + df_flt.norm1
    dep2 = df_flt.support + df_flt.norm2
    df_flt = df_flt[np.logical_and(dep1 >= min_dep, dep2 >= min_dep)]
    print('Filtered out %d SVs based on minimum depth limit' % (n - len(df_flt)))

    n = len(df_flt)
    span = np.array(df_flt.spanning)
    split = np.array(df_flt.split1 + df_flt.split2)
    df_flt = df_flt[np.logical_and(span >= minspan, split >= minsplit)]
    print('Filtered out %d SVs based on spanning/split read limits' % (n - len(df_flt)))

    if filter_chrs:
        n = len(df_flt)
        chr1 = np.array([chrom in valid_chrs for chrom in df_flt.chr1])
        chr2 = np.array([chrom in valid_chrs for chrom in df_flt.chr2])
        df_flt = df_flt[np.logical_and(chr1, chr2)]
        print('Filtered out %d SVs that had non-standard chromosomes' % (n - len(df_flt)))

    if len(blist) > 0:
        n = len(df_flt)
        keep = []
        for idx,sv in df_flt.iterrows():
            pos_olap1 = np.logical_and(sv['pos1'] >= blist.f1, sv['pos1'] <= blist.f2)
            olap1 = np.logical_and(sv['chr1'] == blist.f0, pos_olap1)

            pos_olap2 = np.logical_and(sv['pos2'] >= blist.f1, sv['pos2'] <= blist.f2)
            olap2 = np.logical_and(sv['chr2'] == blist.f0, pos_olap2)

            olaps = blist[np.logical_or(olap1, olap2)]
            if len(olaps) > 0:
                keep.append(False)
            else:
                keep.append(True)

        df_flt = df_flt[keep]
        print('Filtered out %d SVs found in the supplied blacklist' % (n - len(df_flt)))

    return df_flt

def run_simple_snv_filter(snv_df, min_dep, blist, filter_chrs, valid_chrs):

    dep = snv_df['ref'] + snv_df['var']
    snv_df_flt = snv_df[dep>=min_dep]
    print('Filtered out %d SNVs based on minimum depth' % (len(snv_df) - len(snv_df_flt)))

    if filter_chrs:
        chr_flt = np.array([chrom in valid_chrs for chrom in snv_df_flt.chrom])
        snv_df_tmp = snv_df_flt[chr_flt]

        print('Filtered out %d SNVs that had non-standard chromosomes' % (len(snv_df_flt) - len(snv_df_tmp)))
        snv_df_flt = snv_df_tmp

    if len(blist) > 0:
        keep = []
        for idx,snv in snv_df.iterrows():
            pos_olap = np.logical_and(snv['pos']>=blist.f1, snv['pos']<=blist.f2)
            olaps = blist[np.logical_and(snv['chrom']==blist.f0, pos_olap)]
            if len(olaps) > 0:
                keep.append(False)
            else:
                keep.append(True)
        snv_df_flt = snv_df_flt[keep]
        print('Filtered out %d SNVs found in the supplied blacklist' % (len(keep)-sum(keep)))

    return snv_df_flt

def is_clonal_neutral(gtype):
    if gtype=='': return False
    gtype = [x.split(',') for x in gtype.split('|')]
    if len(gtype)==1:
        gt = [float(x) for x in gtype[0]]
        return (gt[0]==1 and gt[1]==1 and gt[2]==1)
    return False

def exceeds_cn_limit(gtype, max_cn):
    '''
    returns True if any one major or minor allele
    of any clone exceeds the max_cn threshold
    '''
    if gtype=='': return False
    gtype = [x.split(',') for x in gtype.split('|')]
    if len(gtype)==1:
        gt = [float(x) for x in gtype[0]]
        return (gt[0] > max_cn or gt[1] > max_cn)
    else:
        return np.any(np.array([float(gt[0]) > max_cn or float(gt[1]) > max_cn for gt in gtype]))

def remove_zero_copynumbers(gtype):
    '''
    remove any clonal or subclonal copy-numbers
    where the total copy-number is zero
    '''
    if gtype=='': return ''
    gtype_tmp = gtype.split('|')
    gtype_tmp = [x.split(',') for x in gtype_tmp]
    if len(gtype)==1:
        gt = [float(x) for x in gtype_tmp[0]]
        if (gt[0]==0 and gt[1]==0) or gt[0] < 0 or gt[1] < 0:
            gtype = ''
    return gtype

def get_weighted_cns(gtypes):
    gtypes_split = [x.split('|') for x in gtypes]
    gtypes_split = [x.split(',') for x in gtypes]
    cn_vals = []
    for gtype in gtypes_split:
        cn_val = sum([int(eval(g[0]))+int(eval(g[1]))*float(eval(g[2])) if g!=[''] else 2 for g in gtype])
        cn_vals.append(cn_val)
    return np.array(cn_vals)/2

def normalise_wins_by_cn(df_flt):
    win1 = df_flt.win_norm1.map(float).values
    win2 = df_flt.win_norm2.map(float).values

    wcn1 = get_weighted_cns(df_flt.gtype1.values)
    wcn2 = get_weighted_cns(df_flt.gtype2.values)

    nonzero1 = np.logical_not(wcn1==0)
    nonzero2 = np.logical_not(wcn2==0)
    win1[nonzero1] = (win1[nonzero1]/wcn1[nonzero1])
    win2[nonzero2] = (win2[nonzero2]/wcn2[nonzero2])

    return win1,win2

def filter_outlying_norm_wins(df_flt):
    win1, win2 = normalise_wins_by_cn(df_flt)
    ranges1 = get_outlier_ranges(win1)
    ranges2 = get_outlier_ranges(win2)

    flt1 = np.logical_and(win1>ranges1[0],win1<ranges1[1])
    flt2 = np.logical_and(win2>ranges2[0],win2<ranges2[1])
    df_flt = df_flt[np.logical_and(flt1,flt2)]

    return df_flt

def run_cnv_filter(df_flt, cnv, ploidy, neutral, filter_outliers, strict_cnv_filt, filter_subclonal, max_cn, are_snvs=False):
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
            is_neutral = [is_clonal_neutral(x) for x in df_flt.gtype.values]
            df_flt2 = df_flt[is_neutral]
            print('Filtered out %d SNVs that were not copy-number neutral' % (n_df - len(df_flt)))

            if filter_outliers:
                n_df = len(df_flt)
                depths = df_flt['ref'].values+df_flt['var'].values
                dep_ranges = get_outlier_ranges(depths)
                df_flt = df_flt[np.logical_and(depths>dep_ranges[0],depths<dep_ranges[1])]
                print('Filtered out %d SNVs which had outlying depths' % (n_df - len(df_flt)))
        else:
            gt1_is_neutral = [is_clonal_neutral(x) for x in df_flt.gtype1.values]
            gt2_is_neutral = [is_clonal_neutral(x) for x in df_flt.gtype2.values]
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
            df_flt['gtype'] = [remove_zero_copynumbers(x) for x in df_flt.gtype.values]

            df_flt = df_flt[df_flt.gtype.values!='']
            print('Filtered out %d SNVs with missing or invalid copy-numbers' % (n_df - len(df_flt)))

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

            if filter_subclonal:
                n_df = len(df_flt)
                is_clon  = [len(x.split('|'))==1 for x in df_flt.gtype.values]
                df_flt = df_flt.loc[is_clon]
                print('Filtered out %d SNVs with subclonal CNV states' % (n_df - len(df_flt)))

            gt_exceeds_cn = np.array([exceeds_cn_limit(gtype, max_cn) for gtype in df_flt.gtype.values])
            if np.any(gt_exceeds_cn):
                n_df = len(df_flt)
                df_flt = df_flt[np.invert(gt_exceeds_cn)]
                print('Filtered out %d SNVs where major or minor alleles exceed the CN max' % (n_df - len(df_flt)))

        else:
            df_flt['gtype1'].fillna('')
            df_flt['gtype2'].fillna('')
            df_flt['gtype1'] = np.array([remove_zero_copynumbers(x) for x in df_flt.gtype1.values])
            df_flt['gtype2'] = np.array([remove_zero_copynumbers(x) for x in df_flt.gtype2.values])

            if strict_cnv_filt:
                # filter out if both CNV states are missing
                df_flt = df_flt[np.logical_and(df_flt.gtype1!='', df_flt.gtype2!='')]
                print('Filtered out %d SVs with missing or invalid copy-numbers' % (n_df - len(df_flt)))
            else:
                # assume normal copy-numbers
                maj_allele = round(ploidy/2) if round(ploidy) > 1 else 1
                min_allele = round(ploidy/2) if round(ploidy) > 1 else 0
                default_gtype = '%d,%d,1.0' % (maj_allele, min_allele)
                df_flt.gtype1[df_flt.gtype1.values==''] = default_gtype
                df_flt.gtype2[df_flt.gtype2.values==''] = default_gtype

            if filter_subclonal:
                n_df = len(df_flt)
                gt1_is_clon  = [len(x.split('|'))==1 for x in df_flt.gtype1.values]
                gt2_is_clon  = [len(x.split('|'))==1 for x in df_flt.gtype2.values]
                is_clonal  = np.logical_and(gt1_is_clon, gt2_is_clon)
                df_flt = df_flt.loc[is_clonal]
                print('Filtered out %d SVs with any subclonal CNV states.' % (n_df - len(df_flt)))

            if filter_outliers:
                n_df = len(df_flt)
                df_flt = filter_outlying_norm_wins(df_flt)
                print('Filtered out %d SVs which had outlying depths' % (n_df - len(df_flt)))

            gt1_exceeds_cn = np.array([exceeds_cn_limit(gtype, max_cn) for gtype in df_flt.gtype1.values])
            gt2_exceeds_cn = np.array([exceeds_cn_limit(gtype, max_cn) for gtype in df_flt.gtype2.values])

            if np.any(gt1_exceeds_cn) or np.any(gt2_exceeds_cn):
                n_df = len(df_flt)
                df_flt = df_flt[np.invert(np.logical_or(gt1_exceeds_cn, gt2_exceeds_cn))]
                print('Filtered out %d SVs where major or minor alleles exceed the CN max' % (n_df - len(df_flt)))

    return df_flt

def match_snv_copy_numbers(snv_df, cnv_df):
    bp_chroms = np.unique(snv_df['chrom'].values)
    bp_chroms = sorted(bp_chroms, key=lambda item: (int(item) if item.isdigit() else str(item)))

    for bchr in bp_chroms:
        gtypes = []
        current_chr = snv_df['chrom'].values==bchr
        var_tmp = snv_df[current_chr]
        cnv_tmp = cnv_df[cnv_df['chr']==bchr]

        if len(cnv_tmp)==0:
            continue

        for pos in var_tmp['pos']:
            cnv_start_list = cnv_tmp.startpos.values
            cnv_end_list   = cnv_tmp.endpos.values
            overlaps = np.logical_and(pos >= cnv_start_list, pos <= cnv_end_list)
            match = cnv_tmp[overlaps]

            if len(match)==0:
                gtypes.append('')
            else:
                gtype = match.loc[match.index[0]].gtype
                gtypes.append(gtype)

        snv_indexes = snv_df[current_chr].index.values
        snv_df.loc[snv_indexes,'gtype'] = gtypes
    return snv_df

def match_copy_numbers(var_df, cnv_df, strict_cnv_filt, sv_offset, bp_fields=['chr1','pos1','dir1','classification','pos2'], gtype_field='gtype1'):

    chrom_field, pos_field, dir_field, class_field, other_pos_field = bp_fields

    var_df[chrom_field]     = [str(x) for x in var_df[chrom_field].values]
    var_df[pos_field]       = [int(x) for x in var_df[pos_field].values]
    var_df[dir_field]       = [str(x) for x in var_df[dir_field].values]
    var_df[class_field]     = [str(x) for x in var_df[class_field].values]
    var_df[other_pos_field] = [int(x) for x in var_df[other_pos_field].values]

    bp_chroms = np.unique(var_df[chrom_field].values)
    bp_chroms = sorted(bp_chroms, key=lambda item: (int(item.partition(' ')[0]) \
                        if item[0].isdigit() else float('inf'), item))

    adj_cnv_field   = '%s_adjacent' % gtype_field
    cnv_dist_field  = '%s_cnv_boundary_dist' % gtype_field

    for bchr in bp_chroms:
        gtypes, cnv_dists, adj_cnvs = [],[],[]
        current_chr = var_df[chrom_field].values==bchr
        var_tmp = var_df[current_chr]
        cnv_tmp = cnv_df[cnv_df['chr']==bchr]

        if len(cnv_tmp)==0:
            var_indexes = var_df[current_chr].index.values
            var_df.loc[var_indexes,gtype_field] = ''
            var_df.loc[var_indexes,adj_cnv_field] = ''
            var_df.loc[var_indexes,cnv_dist_field] = float('nan')
            continue

        sv_info = zip(var_tmp[pos_field],var_tmp[dir_field],var_tmp[class_field],var_tmp[other_pos_field])
        for pos,direct,classification,otherpos in sv_info:
            cnv_gtype,adj_cnv = '',''

            if classification == 'INTRX':
                adjpos = pos - sv_offset if direct == '+' else pos + sv_offset
            else:
                adjpos = pos - sv_offset if gtype_field == 'gtype1' else pos + sv_offset

            adjpos = -1 if direct=='.' else adjpos
            cnv_start_list = cnv_tmp.startpos.values
            cnv_end_list   = cnv_tmp.endpos.values

            closest_start  = min(abs(cnv_start_list-pos))
            closest_end    = min(abs(cnv_end_list-pos))
            cnv_dist       = closest_start if closest_start < closest_end else closest_end

            cnv_dists.append(cnv_dist)

            #print("Checking overlaps for pos %d\n" %adjpos)
            overlaps = np.logical_and(adjpos >= cnv_start_list, adjpos <= cnv_end_list)

            #print(overlaps)
            match = cnv_tmp[overlaps]
            if len(match)==0 and strict_cnv_filt:
                gtypes.append('')
                # add closest CNV to adjacent
                if closest_start < closest_end:
                    start_match = closest_start==abs(cnv_start_list-pos)
                    cnv_gtype   = cnv_tmp[start_match].gtype.values[0]
                else:
                    end_match   = closest_end==abs(cnv_end_list-pos)
                    cnv_gtype   = cnv_tmp[end_match].gtype.values[0]
                adj_cnvs.append(cnv_gtype)
            elif len(match)==0:
                # assign closest CNV state with closest boundary to SNV
                if closest_start < closest_end:
                    start_match = closest_start==abs(cnv_start_list-pos)
                    cnv_gtype   = cnv_tmp[start_match].gtype.values[0]
                    cnv_pos     = cnv_tmp[start_match].startpos.values[0]
                    adj_cnv     = get_adjacent_cnv(cnv_tmp,start_match,pos,cnv_pos,False).gtype.values[0]
                else:
                    end_match   = closest_end==abs(cnv_end_list-pos)
                    cnv_gtype   = cnv_tmp[end_match].gtype.values[0]
                    cnv_pos     = cnv_tmp[end_match].endpos.values[0]
                    adj_cnv     = get_adjacent_cnv(cnv_tmp,end_match,pos,cnv_pos).gtype.values[0]

                gtypes.append(cnv_gtype)
                adj_cnvs.append(adj_cnv)
            else:
                cnv_gtype = match.loc[match.index[0]].gtype
                gtypes.append(cnv_gtype)

                cnv_pos = match.loc[match.index[0]].startpos
                next_cnv = closest_start > closest_end
                if next_cnv:
                    cnv_pos = match.loc[match.index[0]].endpos

                cnv_match = match.index[0]==cnv_tmp.index.values
                adj_cnv = get_adjacent_cnv(cnv_tmp,cnv_match,pos,cnv_pos,next_cnv).gtype.values[0]
                adj_cnvs.append(adj_cnv)

        var_indexes = var_df[current_chr].index.values
        var_df.loc[var_indexes,gtype_field]   = gtypes
        var_df.loc[var_indexes,adj_cnv_field] = adj_cnvs
        var_df.loc[var_indexes,cnv_dist_field] = cnv_dists

    return var_df


def get_adjacent_cnv(cnv,match,pos,cnv_pos,next_cnv=True):

    if not np.any(match):
        return np.empty(0)

    match = match.copy()
    offset = 1 if next_cnv else -1
    where_match = np.where(match)[0][0]

    if where_match+offset<0 or where_match+offset>=len(cnv):
        cnv[match]
    else:
        match[where_match] = False
        match[where_match+offset] = True

    return cnv[match]

def gtypes_match(gtype1,gtype2):
    if gtype1=='' or gtype2=='':
        return False
    gtype1 = [x for x in gtype1.split('|')]
    gtype1 = [float(x) for x in gtype1.split(',')][:2]
    gtype2 = [x for x in gtype2.split('|')]
    gtype2 = [float(x) for x in gtype2.split(',')][:2]
    return np.all(np.array(gtype1)==np.array(gtype2))

def is_same_sv_germline(sv1,sv2,gl_th):
    sv1_chr1, sv1_bp1, sv1_chr2, sv1_bp2 = sv1
    sv2_chr1, sv2_bp1, sv2_chr2, sv2_bp2 = sv2

    if sv1_chr1==sv2_chr1 and sv1_chr2==sv2_chr2:
        if abs(sv1_bp1-sv2_bp1)<gl_th and abs(sv1_bp2-sv2_bp2)<gl_th:
            return True
    if sv1_chr2==sv2_chr1 and sv1_chr1==sv2_chr2:
        if abs(sv1_bp2-sv2_bp1)<gl_th and abs(sv1_bp1-sv2_bp2)<gl_th:
            return True
    return False

def filter_germline(gml_file,sv_df,rlen,insert,gl_th):
    print("Filtering out germline SVs...")
    df_gml = pd.DataFrame(pd.read_csv(gml_file,delimiter='\t',dtype=None,low_memory=False))
    germline = []

    for idx_sv,sv in sv_df.iterrows():
        sv_loc = [str(sv.chr1),int(sv.pos1),str(sv.chr2),int(sv.pos2)]
        same_chrs = np.logical_and(df_gml.chr1.values==sv_loc[0],df_gml.chr2.values==sv_loc[2])
        df_gml_tmp = df_gml[same_chrs]
        for idx_gml,sv_gml in df_gml_tmp.iterrows():
            sv_loc_gml = [str(sv_gml.chr1),int(sv_gml.pos1),str(sv_gml.chr2),int(sv_gml.pos2)]
            if sv_gml.support>0 and is_same_sv_germline(sv_loc,sv_loc_gml,gl_th):
                germline.append(idx_sv)
                break

    print("Filtered out %d SVs that were found in the germline!" % len(germline))
    return sv_df.drop(germline,axis=0)

def adjust_sv_read_counts(sv_df,pi,pl,min_dep,rlen,Config):
    dna_gain_class = Config.get('SVclasses', 'dna_gain_class').split(',')
    dna_loss_class = Config.get('SVclasses', 'dna_loss_class').split(',')
    support_adjust_factor = float(Config.get('FilterParameters', 'support_adjust_factor'))
    filter_subclonal_cnvs = string_to_bool(Config.get('FilterParameters', 'filter_subclonal_cnvs'))

    gt1_sc  = [float(x.split('|')[0].split(',')[2])<1 if x!='' else False for x in sv_df.gtype1.values]
    gt2_sc  = [float(x.split('|')[0].split(',')[2])<1 if x!='' else False for x in sv_df.gtype2.values]
    one_sc  = np.logical_xor(gt1_sc,gt2_sc)

    n = zip(np.array(sv_df.norm1.values),np.array(sv_df.norm2.values))
    s = np.array(sv_df.split1.values+sv_df.split2.values)
    d = np.array(sv_df.spanning.values)
    sup = np.array(d+s,dtype=float)
    Nvar = len(sv_df)
    norm1, norm2 = sv_df.norm1.map(float).values, sv_df.norm2.map(float).values

    try:
        # adjust normal read counts of duplications
        sv_classes = sv_df.classification.values
        dups = np.array([ sv_class in dna_gain_class for idx,sv_class in enumerate(sv_classes) ])

        adjust_factor = 1. - (float(pi) / pl)
        norm1[dups] = [float(n) * adjust_factor for n in norm1[dups]]
        norm2[dups] = [float(n) * adjust_factor for n in norm2[dups]]
    except AttributeError:
        print('Warning, no valid classifications found. SV read counts cannot be adjusted')

    all_indexes = sv_df.index.values
    sup_adjust_factor = 1 + (support_adjust_factor * pi)
    adjusted_support = [int(round(x)) for x in sup * sup_adjust_factor]
    adjusted_depth1 = [int(round(x)) for x in adjusted_support + norm1]
    adjusted_depth2 = [int(round(x)) for x in adjusted_support + norm2]

    adj_norm_means = np.array([np.mean([float(n1), float(n2)]) for n1, n2 in zip(norm1, norm2)])
    sv_df.loc[all_indexes,'adjusted_norm1'] = norm1
    sv_df.loc[all_indexes,'adjusted_norm2'] = norm2
    sv_df.loc[all_indexes,'adjusted_norm_mean'] = adj_norm_means
    sv_df.loc[all_indexes,'adjusted_support'] = adjusted_support
    sv_df.loc[all_indexes,'adjusted_depth1'] = adjusted_depth1
    sv_df.loc[all_indexes,'adjusted_depth2'] = adjusted_depth2
    sv_df.loc[all_indexes,'raw_mean_vaf'] = sup / (sup + sv_df.raw_norm_mean.map(float).values)
    sv_df.loc[all_indexes,'adjusted_vaf1'] = adjusted_support / (norm1 + adjusted_support)
    sv_df.loc[all_indexes,'adjusted_vaf2'] = adjusted_support / (norm2 + adjusted_support)
    sv_df.loc[all_indexes,'adjusted_mean_vaf'] = sup / (adj_norm_means + adjusted_support)

    return sv_df

def string_to_bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def run(args):
    sample      = args.sample
    svs         = args.procd_svs
    gml         = args.germline
    cnvs        = args.cnvs
    out         = args.out
    param_file  = args.param_file
    snvs        = args.snvs
    snv_format  = args.snv_format
    pp_file     = args.pp_file
    cfg         = args.cfg
    blist_file  = args.blist

    Config = configparser.ConfigParser()
    cfg_file = Config.read(cfg)

    if len(cfg_file)==0:
        raise ValueError('No configuration file found')

    max_cn      = int(Config.get('BamParameters', 'max_cn'))
    valid_chrs  = Config.get('ValidationParameters', 'chroms').split(',')
    gl_th       = int(Config.get('FilterParameters', 'germline_threshold'))
    sv_offset   = int(Config.get('FilterParameters', 'sv_offset'))
    minsplit    = int(Config.get('FilterParameters', 'min_split'))
    minspan     = int(Config.get('FilterParameters', 'min_span'))
    sizefilter  = int(Config.get('FilterParameters', 'size_filter'))
    min_dep     = int(Config.get('FilterParameters', 'min_dep'))
    filter_chrs = string_to_bool(Config.get('FilterParameters', 'filter_chroms'))
    neutral     = string_to_bool(Config.get('FilterParameters', 'neutral'))
    filter_otl  = string_to_bool(Config.get('FilterParameters', 'filter_outliers'))
    strict_cnv_filt = string_to_bool(Config.get('FilterParameters', 'strict_cnv_filt'))
    filter_subclonal_cnvs = string_to_bool(Config.get('FilterParameters', 'filter_subclonal_cnvs'))

    out = sample if out == "" else out

    if out!='' and not os.path.exists(out):
        os.makedirs(out)

    pi, ploidy = svp_load.get_purity_ploidy(pp_file, sample, out)
    rlen, insert, insert_std = svp_load.get_read_params(param_file, sample, out)

    blist = pd.DataFrame()
    if blist_file != '' and blist_file.lower() != 'none':
        blist = pd.DataFrame(svp_load.load_blacklist(blist_file))

    if pi < 0 or pi > 1:
        raise ValueError("Tumour purity value not between 0 and 1!")

    sv_df  = pd.DataFrame()
    cnv_df = pd.DataFrame()
    snv_df = pd.DataFrame()

    if snvs!="":
        if snv_format == 'sanger':
            snv_df = load_data.load_snvs_sanger(snvs)
        elif snv_format == 'mutect':
            snv_df = load_data.load_snvs_mutect(snvs,sample)
        elif snv_format == 'mutect_callstats':
            snv_df = load_data.load_snvs_mutect_callstats(snvs)
        elif snv_format == 'consensus':
            snv_df = load_data.load_snvs_consensus(snvs)
        elif snv_format == 'multisnv':
            snv_df = load_data.load_snvs_multisnv(snvs, sample)
        snv_df = run_simple_snv_filter(snv_df, min_dep, blist, filter_chrs, valid_chrs)

    if svs!="":
        sv_df = load_data.load_svs(svs)
        sv_df = run_simple_filter(sv_df,rlen,insert,minsplit,minspan,sizefilter, \
                                  min_dep,filter_chrs,valid_chrs,blist)
        if gml!="":
            sv_df = filter_germline(gml,sv_df,rlen,insert,gl_th)

    if cnvs!="":
        cnv_df = load_data.load_cnvs(cnvs)

        if len(sv_df)>0:
            print('Matching copy-numbers for SVs...')
            sv_df = match_copy_numbers(sv_df,cnv_df,strict_cnv_filt,sv_offset)
            sv_df = match_copy_numbers(sv_df,cnv_df,strict_cnv_filt,sv_offset,\
                    ['chr2','pos2','dir2','classification','pos1'],'gtype2')
            sv_df = run_cnv_filter(sv_df,cnvs,ploidy,neutral,filter_otl,strict_cnv_filt,
                                   filter_subclonal_cnvs,max_cn)

        if len(snv_df)>0:
            print('Matching copy-numbers for SNVs...')
            snv_df = match_snv_copy_numbers(snv_df,cnv_df)
            snv_df = run_cnv_filter(snv_df,cnvs,ploidy,neutral,filter_otl,strict_cnv_filt,
                                    filter_subclonal_cnvs,max_cn,are_snvs=True)
    else:
        print('No CNV input defined, assuming all loci major/minor allele copy-numbers are ploidy/2')
        if len(sv_df)>0:
            maj_allele = round(ploidy/2) if round(ploidy) > 1 else 1
            min_allele = round(ploidy/2) if round(ploidy) > 1 else 0
            default_gtype = '%d,%d,1.0' % (maj_allele, min_allele)
            sv_df['gtype1'] = default_gtype
            sv_df['gtype2'] = default_gtype
            if filter_otl:
                sv_df = run_cnv_filter(sv_df,cnvs,ploidy,neutral,filter_otl,strict_cnv_filt,
                                       filter_subclonal_cnvs,max_cn)
        if len(snv_df)>0:
            maj_allele = round(ploidy/2) if round(ploidy) > 1 else 1
            min_allele = round(ploidy/2) if round(ploidy) > 1 else 0
            default_gtype = '%d,%d,1.0' % (maj_allele, min_allele)
            snv_df['gtype'] = default_gtype
            snv_df = run_cnv_filter(snv_df,cnvs,ploidy,neutral,filter_otl,strict_cnv_filt,
                                    filter_subclonal_cnvs,max_cn,are_snvs=True)

    if len(sv_df)==0 and len(snv_df)==0:
        raise ValueError('No variants found to output!')

    if len(sv_df)>0:
        sv_df.index = range(len(sv_df)) #reindex
        sv_df = adjust_sv_read_counts(sv_df,pi,ploidy,min_dep,rlen,Config)
        sv_df.to_csv('%s/%s_filtered_svs.tsv'%(out,sample),sep='\t',index=False,na_rep='')
        print('Final filtered SV count: %d' % len(sv_df))

    if len(snv_df)>0:
        snv_df.index = range(len(snv_df)) #reindex
        snv_df.to_csv('%s/%s_filtered_snvs.tsv'%(out,sample),sep='\t',index=False,na_rep='')
        print('Final filtered SNV count: %d' % len(snv_df))
