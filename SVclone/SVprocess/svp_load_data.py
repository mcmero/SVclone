import vcf
import numpy as np
import os
from collections import OrderedDict
from . import svp_dtypes as dtypes

def remove_duplicates(svs):
    for idx,row in enumerate(svs):
        #reorder breakpoints based on position or chromosomes
        sv_id, chr1, pos1, dir1, chr2, pos2, dir2, sv_class, oid, opos1, opos2 = row
        if (chr1!=chr2 and chr1>chr2) or (chr1==chr2 and pos1 > pos2):
            svs[idx] = (sv_id, chr2, pos2, dir2, chr1, pos1, dir1, sv_class, oid, opos1, opos2)
    return np.unique(svs)

def load_input_vcf(svin,class_field,use_dir):
    sv_dtype = [s for i,s in enumerate(dtypes.sv_dtype)]

    sv_vcf = vcf.Reader(filename=svin)
    sv_dict = OrderedDict()
    for sv in sv_vcf:

        if sv.FILTER is not None:
            if len(sv.FILTER)>0:
                continue

        sv_dict[sv.ID] = {'CHROM': sv.CHROM, 'POS': sv.POS, 'INFO': sv.INFO,
                          'REF': str(sv.REF[0]), 'ALT': str(sv.ALT[0])}

    svs = np.empty(0, sv_dtype)
    procd = np.empty(0, dtype='<U50')

    for sv_id in sv_dict:
        try:
            sv = sv_dict[sv_id]
            mate_field = 'PARID' if 'PARID' in sv['INFO'] else 'MATEID'
            mate_id = sv['INFO'][mate_field]
            if type(mate_id) == type([]):
                mate_id = mate_id[0]
            mate = sv_dict[mate_id]

            if (sv_id in procd) or (mate_id in procd):
                continue

            dir1, dir2 = '?', '?'
            chr1 = sv['CHROM']
            pos1 = sv['POS']
            chr2 = mate['CHROM']
            pos2 = mate['POS']
            sv_class = sv['INFO'][class_field] if class_field!='' else ''

            if use_dir:
                if sv['ALT'].startswith(']') or sv['ALT'].startswith('['): dir1 = '-'
                if sv['ALT'].endswith('[') or sv['ALT'].endswith(']'): dir1 = '+'
                if mate['ALT'].startswith(']') or mate['ALT'].startswith('['): dir2 = '-'
                if mate['ALT'].endswith('[') or mate['ALT'].endswith(']'): dir2 = '+'

            procd = np.append(procd,[sv_id,mate_id])
            new_id = sv_id.split('_')
            try:
                new_id = int(new_id[0])
            except ValueError:
                new_id = 0
            new_sv = np.array([(new_id, chr1, pos1, dir1, chr2, pos2, dir2, sv_class, sv_id, pos1, pos2)], dtype=sv_dtype)
            svs = np.append(svs,new_sv)
        except KeyError:
            print("SV %s improperly paired or missing attributes"%sv_id)
            continue

    svs['ID'] = range(0,len(svs)) #re-index
    return svs

def load_input_socrates(svin,use_dir,min_mapq,filt_repeats,Config):
    sv_dtype = dtypes.sv_dtype
    pos_field1 = Config.get('SocratesOpts', 'pos1')
    pos_field2 = Config.get('SocratesOpts', 'pos2')
    dir_field1 = Config.get('SocratesOpts', 'dir1')
    dir_field2 = Config.get('SocratesOpts', 'dir2')
    avg_mapq1_field = Config.get('SocratesOpts', 'avg_mapq1')
    avg_mapq2_field = Config.get('SocratesOpts', 'avg_mapq2')
    repeat1_field = Config.get('SocratesOpts', 'repeat1')
    repeat2_field = Config.get('SocratesOpts', 'repeat2')

    soc_in = np.genfromtxt(svin, delimiter='\t', names=True, dtype=None, invalid_raise=False, encoding='utf-8')
    svs = np.empty(0, dtype=sv_dtype)
    filtered_out = 0

    sv_id = 0
    for row in soc_in:
        try:
            bp1 = row[pos_field1].split(':')
            bp2 = row[pos_field2].split(':')
            chr1, pos1 = bp1[0], int(bp1[1])
            chr2, pos2 = bp2[0], int(bp2[1])
            #classification = row['classification']
            if 'normal' in row.dtype.names:
                # has germline info, filter out
                if row['normal']=='normal':
                    continue
            if row[avg_mapq1_field]<min_mapq or row[avg_mapq2_field]<min_mapq:
                filtered_out += 1
                continue
            if filt_repeats!=[]:
                if row[repeat1_field] in filt_repeats and row[repeat2_field] in filt_repeats:
                    filtered_out += 1
                    continue
            add_sv = np.empty(0)

            dir1 = row[dir_field1] if use_dir else '?'
            dir2 = row[dir_field2] if use_dir else '?'

            add_sv = np.array([(sv_id, chr1, pos1, dir1, chr2, pos2, dir2, '', '', pos1, pos2)], dtype=sv_dtype)
            svs = np.append(svs,add_sv)
            sv_id += 1
        except IndexError:
            raise Exception('Supplied Socrates file does not match column names specified in the parameters.py file')

    print('Filtered out %d Socrates SVs, keeping %d SVs' % (filtered_out,len(svs)))
    return remove_duplicates(svs)

def load_input_simple(svin,use_dir,class_field):
    sv_dtype = dtypes.sv_dtype

    sv_tmp = np.genfromtxt(svin, delimiter='\t', names=True, dtype=None, invalid_raise=False, encoding='utf-8')
    sv_tmp = np.atleast_1d(sv_tmp)
    svs = np.empty(0, dtype=sv_dtype)
    sv_id = 0
    for row in sv_tmp:
        chr1 = str(row['chr1'])
        pos1 = int(row['pos1'])
        chr2 = str(row['chr2'])
        pos2 = int(row['pos2'])
        orig_id = str(row['ID']) if 'ID' in sv_tmp.dtype.names else ''
        orig_pos1 = pos1
        orig_pos2 = pos2
        sv_class = row[class_field] if class_field!='' else ''
        add_sv = np.empty(0)

        dir1, dir2 = '?', '?'
        if use_dir:
            dir1 = str(row['dir1'])
            dir2 = str(row['dir2'])

        add_sv = np.array([(sv_id, chr1, pos1, dir1, chr2, pos2, dir2, sv_class, orig_id, orig_pos1, orig_pos2)], dtype=sv_dtype)
        svs = np.append(svs,add_sv)
        sv_id += 1
    return remove_duplicates(svs)

def load_blacklist(blist_file):
    blist = np.genfromtxt(blist_file, delimiter='\t', names=None, dtype=None, invalid_raise=False, encoding='utf-8')
    descr = blist.dtype.descr
    descr[0] = (descr[0][0], '<U10')
    blist = blist.astype(descr)
    if not len(blist) > 0 or not isinstance(blist[0][1], (int, np.integer)) or not isinstance(blist[0][2], (int, np.integer)):
        print('Supplied blacklist is not a valid bed file of intervals')
        return np.empty(0)
    else:
        #remove chr prefixes
        if "chr" in blist[0][0]:
            for idx,row in enumerate(blist):
                blist[idx][0] = blist[idx][0].split('chr')[1]
        return blist

def get_purity_ploidy(pp_file, sample, out):
    '''
    Gets purity/ploidy values from input file,
    if not found, returns defaults. Writes the
    purity/ploidy file to the default loc if it
    doesn't exist.
    '''
    pi      = 1. #default purity
    pl      = 2. #default ploidy

    default_loc = '%s/purity_ploidy.txt' % out
    pp_file = default_loc if pp_file == '' else pp_file

    if os.path.exists(pp_file):
        pur_pl  = np.genfromtxt(pp_file, delimiter='\t', names=True, dtype=None, encoding='utf-8')
        pi      = float(pur_pl['purity'])
        pl      = float(pur_pl['ploidy'])
    else:
        print('WARNING: No purity/ploidy file found. Assuming purity = %f, ploidy = %f' % (pi,pl))

    if pp_file != default_loc:
        with open('%s/purity_ploidy.txt'%out,'w') as outf:
            outf.write("sample\tpurity\tploidy\n")
            outf.write('%s\t%f\t%f\n'%(sample,pi,pl))

    return pi, pl

def get_read_params(params_file, sample, out):
    '''
    Gets read parameter values from input file,
    if not found, returns defaults. Writes the
    purity/ploidy file to the default loc if it
    doesn't exist.
    '''
    rlen    = 100 #default read length
    insert  = 300 #default insert size
    std     = 20

    default_loc =  '%s/read_params.txt' % out
    params_file = default_loc if params_file == '' else params_file

    if os.path.exists(params_file):
        read_params = np.genfromtxt(params_file, delimiter='\t', names=True, dtype=None, invalid_raise=False, encoding='utf-8')
        rlen        = int(read_params['read_len'])
        insert      = float(read_params['insert_mean'])
        std         = float(read_params['insert_std'])
    else:
        print('WARNING: read_params.txt file not found! Assuming read length = %d, mean insert length = %d' % (rlen,insert))

    if params_file != default_loc:
        with open('%s/read_params.txt'%out,'w') as outf:
            outf.write("sample\tread_len\tinsert_mean\tinsert_std\n")
            outf.write('%s\t%f\t%f\t%f\n\n'%(sample,rlen,insert,std))

    return rlen, insert, std
