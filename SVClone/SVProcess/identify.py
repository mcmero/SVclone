from __future__ import print_function
import csv
import ConfigParser
import numpy as np
import vcf
import pysam
import ipdb

from . import load_data
from . import count
from . import svDetectFuncs as svd
from . import bamtools
from .. import dtypes

def classify_event(sv, sv_id, svd_prev_result, prev_sv):
    svd_result = svd.detect(prev_sv, svd_prev_result, sv)
    svd_prev_result, prev_sv = svd_result, sv

    classification = svd.getResultType(svd_result)
    sv_id = sv_id if svd_result[0] == svd.SVtypes.interspersedDuplication else sv_id+1

    return classification, sv_id, svd_prev_result, prev_sv

def has_mixed_evidence(loc_reads, pos, sc_len, threshold):
    split_reads = [count.is_supporting_split_read_lenient(x, pos, threshold*2) for x in loc_reads]
    split_all = loc_reads[np.where(split_reads)[0]]

    if len(split_all) > 0:
        pos, total = sum(split_all['align_start'] < sc_len), len(split_all)
        if pos/float(total) > 0.1 and pos/float(total) < 0.9:
            return True

    return False

def get_dir_span(span):
    is_rev = np.sum(span['is_reverse'])
    assign_dir = '+' if is_rev <= len(span)/2 else '-'
    return assign_dir

def get_dir_split(split, sc_len):
    #align_mean =  np.mean(split['align_start'])
    #assign_dir = '+' if align_mean < sc_len else '-'
    pos, total = sum(split['align_start'] < sc_len), len(split)
    assign_dir = '+' if pos/float(total) > 0.5 else '-'
    return assign_dir

def get_dir(loc_reads, pos, threshold):
    split_reads = [count.is_supporting_split_read_lenient(x, pos, threshold) for x in loc_reads]
    split_all = loc_reads[np.where(split_reads)[0]]
    if len(split_all) > 0:
        dir_split = get_dir_split(split_all, threshold)
        return dir_split
    else:
        return '?'

def get_consensus_align(loc_reads, pos, threshold):
    split_reads = [count.is_supporting_split_read_lenient(x, pos, threshold*2) for x in loc_reads]
    split_all = loc_reads[np.where(split_reads)[0]]

    if len(split_all) != 0:
        ref_starts = [x for x in split_all['ref_start']]
        ref_ends = [x for x in split_all['ref_end']]

        consensus_align_right = max(set(ref_ends), key=ref_ends.count)
        consensus_align_left = max(set(ref_starts), key=ref_starts.count)

        return consensus_align_right, consensus_align_left+1
    else:
        return 0, 0

def get_bp_dir(sv, loc_reads, pos, sc_len, threshold, bp_num):

    bp_dir = 'bp%d_dir' % bp_num
    sv_class = str(sv['classification'])

    ca_right, ca_left = get_consensus_align(loc_reads, pos, threshold)
    ca_right = ca_right if (ca_right-threshold*2 < pos and ca_right+threshold*2 > pos) else 0
    ca_left = ca_left  if (ca_left-threshold*2 < pos  and ca_left+threshold*2 > pos)  else 0

    mixed_right = has_mixed_evidence(loc_reads, ca_right, sc_len, threshold) if ca_right != 0 else False
    mixed_left = has_mixed_evidence(loc_reads, ca_right, sc_len, threshold) if ca_right != 0 else False

    if has_mixed_evidence(loc_reads, pos, sc_len, threshold) or mixed_right or mixed_left:
        sv[bp_dir] = '?'
        sv['classification'] = 'MIXED' if sv_class == '' else sv_class+';MIXED'
    else:
        sv[bp_dir] = get_dir(loc_reads, pos, threshold)
        if sv[bp_dir] == '?':
            sv['classification'] = 'UNKNOWN_DIR' if sv_class == '' else sv_class+';UNKNOWN_DIR'

    return sv, ca_right, ca_left

def retrieve_loc_reads(sv, bam, max_dep, threshold):
    bp_dtype = [('chrom', 'S20'), ('start', int), ('end', int), ('dir', 'S1')]

    sv_id, bp1_chr, bp1_pos, bp1_dir, \
        bp2_chr, bp2_pos, bp2_dir, sv_class = [h[0] for h in dtypes.sv_dtype]
    bp1 = np.array((sv[bp1_chr], sv[bp1_pos]-(threshold*2), \
                    sv[bp1_pos]+(threshold*2), sv[bp1_dir]), dtype=bp_dtype)
    bp2 = np.array((sv[bp2_chr], sv[bp2_pos]-(threshold*2), \
                    sv[bp2_pos]+(threshold*2), sv[bp2_dir]), dtype=bp_dtype)

    bamf = pysam.AlignmentFile(bam, "rb")
    loc1_reads, err_code1 = count.get_loc_reads(bp1, bamf, max_dep)
    loc2_reads, err_code2 = count.get_loc_reads(bp2, bamf, max_dep)
    bamf.close()

    sv_class = str(sv['classification'])
    if err_code1 == 1 or err_code2 == 1:
        sv['classification'] = 'HIDEP' if sv_class == '' else sv_class+';HIDEP'
    elif err_code1 == 2 or err_code2 == 2:
        sv['classification'] = 'READ_FETCH_FAILED' if sv_class == '' \
                                                    else sv_class+';READ_FETCH_FAILED'

    return sv, loc1_reads, loc2_reads, err_code1, err_code2

def get_dir_info(row, bam, max_dep, sc_len, threshold):

    sv, loc1_reads, loc2_reads, err_code1, err_code2 = retrieve_loc_reads(row.copy(), bam, max_dep, threshold)
    if err_code1 != 0 or err_code2 != 0:
        return sv, (0, 0, 0, 0)

    sv, bp1_ca_right, bp1_ca_left = get_bp_dir(sv, loc1_reads, sv['bp1_pos'], sc_len, threshold, 1)
    sv, bp2_ca_right, bp2_ca_left = get_bp_dir(sv, loc2_reads, sv['bp2_pos'], sc_len, threshold, 2)
    consens_aligns = (bp1_ca_right, bp1_ca_left, bp2_ca_right, bp2_ca_left)

    return sv, consens_aligns

def mixed_svs_remain(svs):
    sv_classes = map(lambda x: x.split(';'), svs['classification'])
    sv_classes = np.array([sv for svc in sv_classes for sv in svc])
    return 'MIXED' in sv_classes

def does_break_match(chr1, pos1, chr2, pos2, threshold):
    return chr1 == chr2 and (pos1 + threshold) > pos2 and (pos1 - threshold) < pos2

def get_matching_svs(sv, bp_chr, bp_pos, svs, threshold, mixed=False):
    sv_matches = np.empty(0, dtype=dtypes.sv_dtype)
    which_matches = []
    for row in svs:
        if mixed and row['classification'] != 'MIXED;MIXED':
            continue
        matches = []
        if not np.all(row == sv) and \
            does_break_match(bp_chr, bp_pos, row['bp1_chr'], row['bp1_pos'], threshold):

            sv_matches = np.append(sv_matches, row)
            matches.append(1)
        if not np.all(row == sv) and \
            does_break_match(bp_chr, bp_pos, row['bp2_chr'], row['bp2_pos'], threshold):

            sv_matches = np.append(sv_matches, row)
            matches.append(2)
        if len(matches) > 0:
            which_matches.append(matches)
    return sv_matches, which_matches

def set_dir_class(sv, bp1_dir, bp2_dir, sv_class, new_pos1, new_pos2):
    sv['bp1_dir'] = bp1_dir
    sv['bp2_dir'] = bp2_dir
    sv['classification'] = sv_class
    if new_pos1 != 0:
        sv['bp1_pos'] = new_pos1
    if new_pos2 != 0:
        sv['bp2_pos'] = new_pos2
    return sv

def is_same_sv(sv1, sv2, threshold):
    sv1_chr1, sv1_bp1, sv1_dir1, sv1_chr2, sv1_bp2, sv1_dir2 = sv1
    sv2_chr1, sv2_bp1, sv2_dir1, sv2_chr2, sv2_bp2, sv2_dir2 = sv2

    if sv1_chr1 == sv2_chr1 and sv1_chr2 == sv2_chr2 and \
        sv1_dir1 == sv2_dir1 and sv1_dir2 == sv2_dir2:

        if abs(sv1_bp1-sv2_bp1) < threshold and abs(sv1_bp2-sv2_bp2) < threshold:
            return True
    # same SVs with the order flipped
    if sv1_chr2 == sv2_chr1 and sv1_chr1 == sv2_chr2 and \
        sv1_dir2 == sv2_dir1 and sv1_dir1 == sv2_dir2:

        if abs(sv1_bp2-sv2_bp1) < threshold and abs(sv1_bp1-sv2_bp2) < threshold:
            return True
    return False

def remove_duplicates(svs, threshold):
    to_delete = []
    for idx1, sv1 in enumerate(svs):
        if idx1+1 == len(svs): 
            break
        for idx2, sv2 in enumerate(svs[idx1+1:]):
            if idx1 != idx2:
                sv1_tmp = (sv1['bp1_chr'], sv1['bp1_pos'], sv1['bp1_dir'],\
                           sv1['bp2_chr'], sv1['bp2_pos'], sv1['bp2_dir'])
                sv2_tmp = (sv2['bp1_chr'], sv2['bp1_pos'], sv2['bp1_dir'],\
                           sv2['bp2_chr'], sv2['bp2_pos'], sv1['bp2_dir'])
                if is_same_sv(sv1_tmp, sv2_tmp, threshold):
                    to_delete.append(idx1)
    svs = np.delete(svs, np.array(to_delete))
    return svs

def get_sv_pos_ranks(sv_list, threshold):
    '''
    return ranks for the positions of a list of svs (in order of pos and sv)
    returns an empty list if the svs are not all on the same chromosome
    '''
    chr_list = [[x['bp1_chr'], x['bp2_chr']] for x in sv_list]
    pos_list = [[x['bp1_pos'], x['bp2_pos']] for x in sv_list]

    ranks = []
    if len(set([c for chrom in chr_list for c in chrom])) == 1:

        all_pos = [p for pos in pos_list for p in pos]
        uniq_pos = []
        for idx, pos in enumerate(all_pos):
            if not np.any(abs(all_pos[idx+1:]-pos) <= threshold):
                uniq_pos.append(pos)
        uniq_pos = np.sort(np.array(uniq_pos))

        for pos in pos_list:
            rank1 = int(np.where(min(abs(uniq_pos-pos[0])) == abs(uniq_pos-pos[0]))[0]+1)
            rank2 = int(np.where(min(abs(uniq_pos-pos[1])) == abs(uniq_pos-pos[1]))[0]+1)
            ranks.append([rank1, rank2])

    return ranks

def split_dirs_dual_mixed_sv(svs, row, idx, ca, threshold):

    bp1_svmatch, bp1_which = get_matching_svs(svs[idx], svs[idx]['bp1_chr'], \
                                svs[idx]['bp1_pos'], svs, threshold, mixed=True)
    bp2_svmatch, bp2_which = get_matching_svs(svs[idx], svs[idx]['bp2_chr'], \
                                svs[idx]['bp2_pos'], svs, threshold, mixed=True)

    if len(bp1_svmatch) > 0 and len(bp2_svmatch) > 0:

        other_idx1 = int(bp1_svmatch[0]['ID'])
        other_idx2 = int(bp2_svmatch[0]['ID'])

        if len(bp1_svmatch) > 1 or len(bp2_svmatch) > 1:
            indexes = [int(sv_match['ID']) for sv_match in bp1_svmatch] + [row['ID']]
            indexes = [int(sv_match['ID']) for sv_match in bp2_svmatch] + indexes

            for tmp_idx in indexes:
                svs[tmp_idx]['classification'] = 'COMPLEX'

            return svs

        if bp1_svmatch[0] == bp2_svmatch[0]:
            # likely an inversion - both sides match same break
            # set one SV to + dirs, the other to - dirs

            new_sv = svs[idx].copy()
            new_pos1, new_pos2 = ca[idx]['bp1_ca_right'], ca[idx]['bp2_ca_right']
            svs[idx] = set_dir_class(new_sv, '+', '+', '', new_pos1, new_pos2)

            new_pos1, new_pos2 = ca[idx]['bp1_ca_left'], ca[idx]['bp2_ca_left']
            svs[other_idx1] = set_dir_class(svs[idx], '-', '-', '', new_pos1, new_pos2)

        elif len(bp1_which[0]) == 1 and len(bp2_which[0]) == 1:
            # break partners differ, 1 side matched - check for translocation

            sv1, sv2, sv3 = row, bp1_svmatch[0], bp2_svmatch[0]
            sv_list = [sv1, sv2, sv3]

            def contains_all_ranks(ranks_to_test, ranks):
                contains_ranks = [np.any(np.all(np.array(ranks) == rank, axis=1)) \
                                    for rank in ranks_to_test]
                return np.all(contains_ranks)

            ranks = get_sv_pos_ranks(sv_list, threshold)
            indexes = [idx, other_idx1, other_idx2]
            transloc_ranks = [[1, 2], [2, 3], [1, 3]]

            if len(ranks) > 0 and contains_all_ranks(transloc_ranks, ranks):
                # translocation signature - adjust points
                for rank, sv, idxs in zip(ranks, sv_list, indexes):
                    dirs = ['-', '+'] if rank == [1, 3] else ['+', '-']
                    new_pos1 = ca[idxs]['bp1_ca_right'] if \
                                dirs[0] == '+' else ca[idxs]['bp1_ca_left']
                    new_pos2 = ca[idxs]['bp2_ca_right'] if \
                                dirs[1] == '+' else ca[idxs]['bp2_ca_left']
                    svs[idxs] = set_dir_class(sv, dirs[0], dirs[1], '', new_pos1, new_pos2)

    elif len(bp1_svmatch) == 0 and len(bp2_svmatch) == 0:
        # probably an inversion that has no partner - not called by SV caller?
        new_pos1, new_pos2 = ca[idx]['bp1_ca_left'], ca[idx]['bp2_ca_left']
        svs[idx] = set_dir_class(svs[idx], '-', '-', '', new_pos1, new_pos2)

        # create a new entry for the inversion partner
        new_sv = svs[idx].copy()
        new_pos1, new_pos2 = ca[idx]['bp1_ca_right'], ca[idx]['bp2_ca_right']
        new_sv = set_dir_class(new_sv, '+', '+', '', new_pos1, new_pos2)
        new_sv['ID'] = svs[len(svs)-1]['ID']+1
        svs = np.append(svs, new_sv)

    else:
        # multiple matching events or only one matched - don't know how to deal with this
        svs[idx]['classification'] = 'UNKNOWN_DIR'

    return svs

def split_mixed_svs(svs, ca, threshold):
    for idx, row in enumerate(svs):
        sv_class = np.array(row['classification'].split(';'))

        if 'MIXED' in sv_class:

            if 'UNKNOWN_DIR' in sv_class:
                svs[idx]['classification'] = 'UNKNOWN_DIR'
                break #can't fix this if one direction is unknown

            if len(sv_class) > 1 and sv_class[0] == 'MIXED' and sv_class[1] == 'MIXED':
                svs = split_dirs_dual_mixed_sv(svs, row, idx, ca, threshold)

            elif row['bp1_dir'] != '?':
                #set dir to -, create new sv with dir set to +
                new_sv = svs[idx].copy()
                new_pos2 = ca[idx]['bp2_ca_right']
                new_sv = set_dir_class(new_sv, row['bp1_dir'], '+', '', 0, new_pos2)
                svs = np.append(svs, new_sv)

                new_pos2 = ca[idx]['bp2_ca_left']
                svs[idx] = set_dir_class(svs[idx], row['bp1_dir'], '-', '', 0, new_pos2)

            elif row['bp2_dir'] != '?':
                #set dir to -, create new sv with dir set to +
                new_sv = svs[idx].copy()
                new_pos1 = ca[idx]['bp1_ca_right']
                new_sv = set_dir_class(new_sv, '+', row['bp2_dir'], '', new_pos1, 0)
                svs = np.append(svs, new_sv)

                svs[idx] = set_dir_class(svs[idx], '-', row['bp2_dir'], '', new_pos1, 0)

            else:
                svs[idx]['classification'] = 'UNKNOWN_DIR'

            break #can only fix one set at a time

    return svs

def preproc_svs(args):

    svin         = args.svin
    bam          = args.bam
    out          = args.out
    simple       = args.simple_svs
    socrates     = args.socrates
    use_dir      = args.use_dir
    min_mapq     = args.min_mapq
    class_field  = args.class_field
    filt_repeats = args.filt_repeats
    trust_sc_pos = args.trust_sc_pos
    rlen         = args.rlen

    cfg = args.cfg
    Config = ConfigParser.ConfigParser()
    cfg_file = Config.read(cfg)

    if len(cfg_file)==0:
        raise ValueError('No configuration file found')

    max_cn   = int(Config.get('GlobalParameters', 'max_cn'))
    mean_cov = int(Config.get('GlobalParameters', 'mean_cov'))
    sc_len   = int(Config.get('GlobalParameters', 'sc_len'))
    threshold = int(Config.get('GlobalParameters', 'threshold'))

    filt_repeats = filt_repeats.split(', ') if filt_repeats != '' else filt_repeats
    filt_repeats = [rep for rep in filt_repeats if rep != '']

    outname = '%s_svin.txt'%out

    svs = np.empty(0)
    if simple:
        svs = load_data.load_input_simple(svin, use_dir, class_field)
    elif socrates:
        svs = load_data.load_input_socrates(svin, use_dir, min_mapq, filt_repeats, Config)
    else:
        svs = load_data.load_input_vcf(svin, class_field)

    consens_dtype = [('bp1_ca_right', int), ('bp1_ca_left', int), \
                        ('bp2_ca_right', int), ('bp2_ca_left', int)]
    ca = np.zeros(len(svs), dtype=consens_dtype)

    if rlen<0:
        rlen = bamtools.estimateTagSize(bam)
    
    inserts = bamtools.estimateInsertSizeDistribution(bam)
    inserts = (max(rlen*2,inserts[0]),inserts[1])
    
    max_ins = inserts[0]+(3*inserts[1]) #max fragment size = mean fragment len + (fragment std * 3)
    max_dep = ((mean_cov*(max_ins*2))/rlen)*max_cn

    if not use_dir:
        for idx, sv in enumerate(svs):
            svs[idx], ca[idx] = get_dir_info(sv, bam, max_dep, sc_len, threshold)

        while mixed_svs_remain(svs):
            svs = split_mixed_svs(svs, ca, threshold)
            print('%d dual mixed classifications remain' % \
                sum(svs['classification'] == 'MIXED;MIXED'))

    elif not trust_sc_pos:
        # use directions, but recheck the consensus alignments
        for idx, sv in enumerate(svs):
            sv_tmp, loc1_reads, loc2_reads, err_code1, err_code2 = \
                retrieve_loc_reads(sv.copy(), bam, max_dep, threshold)

            if err_code1 != 0 or err_code2 != 0:
                continue

            ca_right, ca_left = get_consensus_align(loc1_reads, sv['bp1_pos'], threshold)
            new_align = ca_right if sv['bp1_dir'] == '+' else ca_left
            if new_align != 0:
                svs[idx]['bp1_pos'] = new_align

            ca_right, ca_left = get_consensus_align(loc2_reads, sv['bp2_pos'], threshold)
            new_align = ca_right if sv['bp2_dir'] == '+' else ca_left
            if new_align != 0:
                svs[idx]['bp2_pos'] = new_align

    print('Removing duplicate SVs...')
    svs = remove_duplicates(svs, threshold)

    sv_id = 0
    svd_prev_result, prev_sv = None, None
    for idx, sv in enumerate(svs):
        if sv['classification'] == '':
            sv_class, sv_id, svd_prev_result, prev_sv = \
                classify_event(sv, sv_id, svd_prev_result, prev_sv)
            svs[idx]['classification'] = sv_class
        else:
            sv_id += 1
        svs[idx]['ID'] = sv_id

    # reclassify once all have been processed
    for idx, sv in enumerate(svs):
        trx_label = svd.getResultType([svd.SVtypes.translocation])
        intins_label = svd.getResultType([svd.SVtypes.interspersedDuplication])
        if sv['classification'] == intins_label:
            svs[idx-1]['classification'] = intins_label
            translocs = svd.detectTransloc(idx, svs, threshold)
            if len(translocs) > 0:
                new_idx = idx-1
                for i in translocs:
                    svs[i]['classification'] = trx_label
                    svs[i]['ID'] = new_idx

    with open(outname, 'w') as outf:
        writer = csv.writer(outf, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
        header = [field for field, dtype in dtypes.sv_dtype]
        writer.writerow(header)
        for sv_out in svs:
            writer.writerow(sv_out)
