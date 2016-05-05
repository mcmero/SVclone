from __future__ import print_function
import csv
import ConfigParser
import numpy as np
import vcf
import pysam
import os

from . import load_data
from . import count
from . import svDetectFuncs as svd
from . import bamtools
from . import dtypes

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

    bp_dir = 'dir%d' % bp_num
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

    sv_id, chr1, pos1, dir1, \
        chr2, pos2, dir2, sv_class = [h[0] for h in dtypes.sv_dtype]
    bp1 = np.array((sv[chr1], sv[pos1]-(threshold*2), \
                    sv[pos1]+(threshold*2), sv[dir1]), dtype=bp_dtype)
    bp2 = np.array((sv[chr2], sv[pos2]-(threshold*2), \
                    sv[pos2]+(threshold*2), sv[dir2]), dtype=bp_dtype)

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

    sv, ca_right1, ca_left1 = get_bp_dir(sv, loc1_reads, sv['pos1'], sc_len, threshold, 1)
    sv, ca_right2, ca_left2 = get_bp_dir(sv, loc2_reads, sv['pos2'], sc_len, threshold, 2)
    consens_aligns = (ca_right1, ca_left1, ca_right2, ca_left2)

    return sv, consens_aligns

def num_mixed_svs(svs):
    sv_classes = map(lambda x: x.split(';'), svs['classification'])
    sv_classes = np.array([sv for svc in sv_classes for sv in svc])
    sv_classes = sv_classes = np.array(['MIXED' in svc for svc in sv_classes])
    return sum(sv_classes)

def does_break_match(chr1, pos1, chr2, pos2, threshold):
    return chr1 == chr2 and (pos1 + threshold) > pos2 and (pos1 - threshold) < pos2

def get_matching_svs(idx, sv, bp_chr, bp_pos, svs, threshold, mixed=False):
    sv_matches = np.empty(0, dtype=dtypes.sv_dtype)
    which_matches = []

    for i,row in enumerate(svs):
        if mixed and row['classification'] != 'MIXED;MIXED':
            continue
        if i==idx:
            continue
        matches = []

        if does_break_match(bp_chr, bp_pos, row['chr1'], row['pos1'], threshold):
            match_sv = np.array(row, dtype=dtypes.sv_dtype)
            sv_matches = np.append(sv_matches, match_sv)
            matches.append(1)

        if does_break_match(bp_chr, bp_pos, row['chr2'], row['pos2'], threshold):
            match_sv = np.array(row, dtype=dtypes.sv_dtype)
            sv_matches = np.append(sv_matches, match_sv)
            matches.append(2)

        if len(matches) > 0:
            which_matches.append(matches)

    return sv_matches, which_matches

def set_dir_class(sv, dir1, dir2, sv_class, new_pos1, new_pos2):
    sv['dir1'] = dir1
    sv['dir2'] = dir2
    sv['classification'] = sv_class
    if new_pos1 != 0:
        sv['pos1'] = new_pos1
    if new_pos2 != 0:
        sv['pos2'] = new_pos2
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
                sv1_tmp = (sv1['chr1'], sv1['pos1'], sv1['dir1'],\
                           sv1['chr2'], sv1['pos2'], sv1['dir2'])
                sv2_tmp = (sv2['chr1'], sv2['pos1'], sv2['dir1'],\
                           sv2['chr2'], sv2['pos2'], sv1['dir2'])
                if is_same_sv(sv1_tmp, sv2_tmp, threshold):
                    to_delete.append(idx1)
    svs = np.delete(svs, np.array(to_delete))
    svs['ID'] = range(0,len(svs)) #re-index
    return svs

def get_sv_pos_ranks(sv_list, threshold):
    '''
    return ranks for the positions of a list of svs (in order of pos and sv)
    returns an empty list if the svs are not all on the same chromosome
    '''
    chr_list = [[x['chr1'], x['chr2']] for x in sv_list]
    pos_list = [[x['pos1'], x['pos2']] for x in sv_list]

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

def set_svs_as_complex(svs, svmatch):
    # ensure that svs are of the correct dtype
    svs = np.array(svs, dtype=dtypes.sv_dtype)
    svmatch = np.array(svmatch, dtype=dtypes.sv_dtype)
    for match in svmatch:
        for tmp_id in np.where(svs==match)[0]:
            svs[tmp_id]['classification'] = 'COMPLEX'
    return svs

def split_dirs_dual_mixed_sv(svs, row, idx, ca, threshold):

    svmatch1, which1 = get_matching_svs(idx, svs[idx], svs[idx]['chr1'], \
                                svs[idx]['pos1'], svs, threshold, mixed=True)
    svmatch2, which2 = get_matching_svs(idx, svs[idx], svs[idx]['chr2'], \
                                svs[idx]['pos2'], svs, threshold, mixed=True)

    if len(svmatch1) > 0 and len(svmatch2) > 0:
        
        if len(svmatch1) > 1 or len(svmatch2) > 1:

            svs = set_svs_as_complex(svs, svmatch1)
            svs = set_svs_as_complex(svs, svmatch2)
            return svs

        # ensure that svs are of the correct dtype
        svs = np.array(svs, dtype=dtypes.sv_dtype)
        other_idx1 = np.where(svs==svmatch1)[0][0]
        other_idx2 = np.where(svs==svmatch2)[0][0]

        if svmatch1[0] == svmatch2[0]:
            # likely an inversion - both sides match same break
            # set one SV to + dirs, the other to - dirs

            new_sv = svs[idx].copy()
            new_pos1, new_pos2 = ca[idx]['ca_right1'], ca[idx]['ca_right2']
            svs[idx] = set_dir_class(new_sv, '+', '+', '', new_pos1, new_pos2)

            new_pos1, new_pos2 = ca[idx]['ca_left1'], ca[idx]['ca_left2']
            svs[other_idx1] = set_dir_class(svs[idx], '-', '-', '', new_pos1, new_pos2)

        elif len(which1[0]) == 1 and len(which2[0]) == 1:
            # break partners differ, 1 side matched - check for translocation

            sv1, sv2, sv3 = row, svmatch1[0], svmatch2[0]
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
                    new_pos1 = ca[idxs]['ca_right1'] if \
                                dirs[0] == '+' else ca[idxs]['ca_left1']
                    new_pos2 = ca[idxs]['ca_right2'] if \
                                dirs[1] == '+' else ca[idxs]['ca_left2']
                    svs[idxs] = set_dir_class(sv, dirs[0], dirs[1], '', new_pos1, new_pos2)
            else:                
                svs = set_svs_as_complex(svs, svmatch1)
                svs = set_svs_as_complex(svs, svmatch2)
                return svs

    elif len(svmatch1) == 0 and len(svmatch2) == 0:
        # probably an inversion that has no partner - not called by SV caller?
        new_pos1, new_pos2 = ca[idx]['ca_left1'], ca[idx]['ca_left2']
        svs[idx] = set_dir_class(svs[idx], '-', '-', '', new_pos1, new_pos2)

        # create a new entry for the inversion partner
        new_sv = svs[idx].copy()
        new_pos1, new_pos2 = ca[idx]['ca_right1'], ca[idx]['ca_right2']
        new_sv = set_dir_class(new_sv, '+', '+', '', new_pos1, new_pos2)
        new_sv['ID'] = svs[len(svs)-1]['ID']+1
        svs = np.append(svs, new_sv)

    else:
        # multiple matching events or only one matched - don't know how to deal with this
        svs[idx]['classification'] = 'UNKNOWN_DIR'

    return svs

def split_mixed_svs(svs, ca, threshold):
    new_svs = []
    for idx, row in enumerate(svs):

        sv_class = np.array(row['classification'].split(';'))
        if 'MIXED' in sv_class:

            if 'UNKNOWN_DIR' in sv_class:
                svs[idx]['classification'] = 'UNKNOWN_DIR'
                break #can't fix this if one direction is unknown

            elif len(sv_class) > 1 and sv_class[0] == 'MIXED' and sv_class[1] == 'MIXED':
                new_svs = split_dirs_dual_mixed_sv(svs, row, idx, ca, threshold)

            elif row['dir1'] in ['-','+']:
                #set dir to -, create new sv with dir set to +                
                new_sv = svs[idx].copy()
                new_pos2 = ca[idx]['ca_right2']
                new_sv = set_dir_class(new_sv, row['dir1'], '+', '', 0, new_pos2)
                new_sv['ID'] = svs[len(svs)-1]['ID']+1
                svs[idx] = set_dir_class(svs[idx], row['dir1'], '-', '', 0, new_pos2)

                new_svs = np.append(svs, new_sv)

            elif row['dir2'] in ['-','+']:
                #set dir to -, create new sv with dir set to +
                new_sv = svs[idx].copy()
                new_pos1 = ca[idx]['ca_right1']
                new_sv = set_dir_class(new_sv, '+', row['dir2'], '', new_pos1, 0)
                new_sv['ID'] = svs[len(svs)-1]['ID']+1
                svs[idx] = set_dir_class(svs[idx], '-', row['dir2'], '', new_pos1, 0)

                new_svs = np.append(svs, new_sv)

            else:
                svs[idx]['classification'] = 'UNKNOWN_DIR'

            break #can only fix one set at a time

    if len(new_svs)>0:
        return new_svs
    else:
        return svs

def sv_in_blacklist(sv, blist):

    for row in blist:
        if sv['chr1'] == row[0]:
            if sv['pos1'] >= row[1] and row[2] >= sv['pos1']:
                return True
        if sv['chr2'] == row[0]:
            if sv['pos2'] >= row[1] and row[2] >= sv['pos2']:
                return True

    return False

def infer_sv_dirs(svs, ca, bam, max_dep, sc_len, threshold, blist):

    print('Inferring SV directions...')
    for idx, sv in enumerate(svs):
        if len(blist) > 0 and sv_in_blacklist(sv, blist):
            svs[idx]['classification'] = 'BLACKLIST'
            continue
        svs[idx], ca[idx] = get_dir_info(sv, bam, max_dep, sc_len, threshold)

#        tmp_out = '%s_dirout.txt' % out
#        svs = np.genfromtxt(tmp_out, delimiter='\t', names=True, dtype=None, invalid_raise=False)
         # write directionality output (for debugging)
#        with open(tmp_out, 'w') as outf:
#            writer = csv.writer(outf, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar='')
#            header = [field for field, dtype in dtypes.sv_dtype]
#            writer.writerow(header)
#            for sv_out in svs:
#                writer.writerow(sv_out)

    while num_mixed_svs(svs)>0:

        before = num_mixed_svs(svs)
        svs = split_mixed_svs(svs, ca, threshold)
        after = num_mixed_svs(svs)
        print('%d mixed classifications remain' % after)

        # to prevent infinite loops...
        if before == after:
            print('Failed to fix mixed direction SVs. Skipping...')
            break

    return svs, ca

def string_to_bool(v):
  return v.lower() in ("yes", "true", "t", "1")

def preproc_svs(args):

    svin         = args.svin
    bam          = args.bam
    out          = args.out
    sample       = args.sample
    sv_format    = args.sv_format
    blist_file   = args.blist

    cfg = args.cfg
    Config = ConfigParser.ConfigParser()
    cfg_file = Config.read(cfg)

    if len(cfg_file)==0:
        raise ValueError('No configuration file found')

    max_cn       = int(Config.get('BamParameters', 'max_cn'))
    mean_cov     = int(Config.get('BamParameters', 'mean_cov'))
    rlen         = int(Config.get('BamParameters', 'read_len'))

    use_dir      = string_to_bool(Config.get('SVidentifyParameters', 'use_dir'))
    trust_sc_pos = string_to_bool(Config.get('SVidentifyParameters', 'trust_sc_position'))
    sc_len       = int(Config.get('SVcountParameters', 'sc_len'))
    threshold    = int(Config.get('SVcountParameters', 'threshold'))

    min_mapq     = int(Config.get('SocratesOpts', 'min_mapq'))
    filt_repeats = str(Config.get('SocratesOpts', 'filter_repeats'))
    filt_repeats = '' if filt_repeats == 'none' else filt_repeats
    filt_repeats = filt_repeats.split(', ') if filt_repeats != '' else filt_repeats
    filt_repeats = [rep for rep in filt_repeats if rep != '']

    class_field  = str(Config.get('SVidentifyParameters', 'sv_class_field'))
    class_field  = '' if class_field == 'none' else class_field

    out = sample if out == "" else out
    outname = '%s/%s_svin.txt' % (out, sample)

    if out!='' and not os.path.exists(out):
        os.makedirs(out)

    svs = np.empty(0)
    print('Loading SV calls...')
    if sv_format == 'vcf':
        svs = load_data.load_input_vcf(svin, class_field, use_dir)
    elif sv_format == 'simple':
        svs = load_data.load_input_simple(svin, use_dir, class_field)
    elif sv_format == 'socrates':
        svs = load_data.load_input_socrates(svin, use_dir, min_mapq, filt_repeats, Config)
    else:
        raise ValueError('Valid input format not specified.')

    blist = np.empty(0)
    if blist_file != '':
        blist = load_data.load_blacklist(blist_file)

    consens_dtype = [('ca_right1', int), ('ca_left1', int), \
                        ('ca_right2', int), ('ca_left2', int)]
    ca = np.zeros(len(svs), dtype=consens_dtype)

    if rlen<0:
        rlen = bamtools.estimateTagSize(bam)
    
    inserts = bamtools.estimateInsertSizeDistribution(bam)
    inserts = (max(rlen*2,inserts[0]),inserts[1])
    
    max_ins = inserts[0]+(3*inserts[1]) #max fragment size = mean fragment len + (fragment std * 3)
    max_dep = ((mean_cov*(max_ins*2))/rlen)*max_cn
    
    if not use_dir:
        svs, ca = infer_sv_dirs(svs, ca, bam, max_dep, sc_len, threshold, blist)
    elif not trust_sc_pos:
        print('Recalibrating consensus alignments...')
        # use directions, but recheck the consensus alignments
        for idx, sv in enumerate(svs):
            if len(blist) > 0 and sv_in_blacklist(sv, blist):
                svs[idx]['classification'] = 'BLACKLIST'
                continue

            sv_tmp, loc1_reads, loc2_reads, err_code1, err_code2 = \
                retrieve_loc_reads(sv.copy(), bam, max_dep, threshold)

            if err_code1 != 0 or err_code2 != 0:
                continue

            ca_right, ca_left = get_consensus_align(loc1_reads, sv['pos1'], threshold)
            new_align = ca_right if sv['dir1'] == '+' else ca_left
            if new_align != 0:
                svs[idx]['pos1'] = new_align

            ca_right, ca_left = get_consensus_align(loc2_reads, sv['pos2'], threshold)
            new_align = ca_right if sv['dir2'] == '+' else ca_left
            if new_align != 0:
                svs[idx]['pos2'] = new_align

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
