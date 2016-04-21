import numpy as np

#####################################################################
# Data dtypes
#####################################################################

bp_dtype = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S1')]

sv_dtype = [('ID', 'int64'),
            ('bp1_chr', 'S20'),
            ('bp1_pos', 'int64'),
            ('bp1_dir', 'S1'),
            ('bp2_chr', 'S20'),
            ('bp2_pos', 'int64'),
            ('bp2_dir', 'S1'),
            ('classification', 'S100')] 

read_dtype = [('query_name', 'S150'), 
            ('chrom', 'S50'), 
            ('ref_start', 'int64'), 
            ('ref_end', 'int64'), 
            ('align_start', 'int64'), 
            ('align_end', 'int64'), 
            ('len', 'int64'), 
            ('ins_len', 'int64'), 
            ('is_reverse', np.bool)]

#####################################################################
# Output dtypes
#####################################################################

# Output from count step
sv_out_dtype = [('ID', 'int64'),
            ('bp1_chr', 'S20'),
            ('bp1_pos', 'int64'),
            ('bp1_dir', 'S1'),
            ('bp2_chr', 'S20'),
            ('bp2_pos', 'int64'),
            ('bp2_dir', 'S1'),
            ('classification', 'S100'),
            ('bp1_split_norm', 'int64'),
            ('bp1_norm_olap_bp', 'int64'),
            ('bp1_span_norm', 'int64'),
            ('bp1_win_norm', 'int64'), 
            ('bp1_split', 'int64'), 
            ('bp1_sc_bases', 'int64'), 
            ('bp1_total_reads', 'int64'),
            ('bp2_split_norm', 'int64'), 
            ('bp2_norm_olap_bp', 'int64'),
            ('bp2_span_norm', 'int64'), 
            ('bp2_win_norm', 'int64'), 
            ('bp2_split', 'int64'), 
            ('bp2_sc_bases', 'int64'), 
            ('bp2_total_reads', 'int64'),
            ('anomalous', 'int64'),
            ('spanning', 'int64'),
            ('norm1', 'int64'),
            ('norm2', 'int64'), 
            ('support', 'int64'),
            ('vaf1', float),
            ('vaf2', float)]

# SV copy-number output
sv_cn_dtype = [('chr1', 'S50'),
            ('pos1', int),
            ('total_copynumber1', int),
            ('no_chrs_bearing_mutation1', int),
            ('chr2', 'S50'),
            ('pos2', int),
            ('total_copynumber2', int),
            ('no_chrs_bearing_mutation2', int)]

# SV most likely copy-numbers output
sv_mlcn_dtype =[('chr1', 'S50'),
            ('pos1', int),
            ('chr2', 'S50'),
            ('pos2', int),
            ('pos1_bb_CN', 'S50'),
            ('pos2_bb_CN', 'S50'),
            ('most_likely_ref_copynumber', int),
            ('most_likely_variant_copynumber', int),
            ('prop_chrs_bearing_mutation', float),
            ('support', int),
            ('depth', int),      
            ('pv', float)]

# SV multiplicities output
sv_mult_dtype =[('chr1', 'S50'),
            ('pos1', int),
            ('chr2', 'S50'),
            ('pos2', int),
            ('tumour_copynumber', int),
            ('multiplicity', int),
            ('tumour_copynumber_options', 'S50'),
            ('multiplicity_options', 'S50'),
            ('probabilities', 'S50')]

# SNV copy-number output
snv_cn_dtype =  [('chr', 'S50'),
            ('pos', int),
            ('total_copynumber', int),
            ('no_chrs_bearing_mutation', int)]

# SNV most likely copy-number output
snv_mlcn_dtype =[('chr', 'S50'),
            ('pos', int),
            ('bb_CN', 'S50'),
            ('most_likely_ref_copynumber', int),
            ('most_likely_variant_copynumber', int),
            ('prop_chrs_bearing_mutation', float)]

# SNV multiplicities output
snv_mult_dtype =[('chr', 'S50'),
            ('pos', int),
            ('tumour_copynumber', int),
            ('multiplicity', int),
            ('tumour_copynumber_options', 'S50'),
            ('multiplicity_options', 'S50'),
            ('probabilities', 'S50')]
