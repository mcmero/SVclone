import numpy as np

#####################################################################
# Data dtypes
#####################################################################

bp_dtype = [('chrom','<U20'),('start', int), ('end', int), ('dir', '<U1')]

sv_dtype = [('ID', 'int64'),
            ('chr1', '<U20'),
            ('pos1', 'int64'),
            ('dir1', '<U1'),
            ('chr2', '<U20'),
            ('pos2', 'int64'),
            ('dir2', '<U1'),
            ('classification', '<U100')] 

read_dtype = [('query_name', '<U150'), 
            ('chrom', '<U50'), 
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
            ('chr1', '<U20'),
            ('pos1', 'int64'),
            ('dir1', '<U1'),
            ('chr2', '<U20'),
            ('pos2', 'int64'),
            ('dir2', '<U1'),
            ('classification', '<U100'),
            ('split_norm1', 'int64'),
            ('norm_olap_bp1', 'int64'),
            ('span_norm1', 'int64'),
            ('win_norm1', 'int64'), 
            ('split1', 'int64'), 
            ('sc_bases1', 'int64'), 
            ('total_reads1', 'int64'),
            ('split_norm2', 'int64'), 
            ('norm_olap_bp2', 'int64'),
            ('span_norm2', 'int64'), 
            ('win_norm2', 'int64'), 
            ('split2', 'int64'), 
            ('sc_bases2', 'int64'), 
            ('total_reads2', 'int64'),
            ('anomalous', 'int64'),
            ('spanning', 'int64'),
            ('norm1', 'int64'),
            ('norm2', 'int64'), 
            ('support', 'int64'),
            ('vaf1', float),
            ('vaf2', float)]

# SV copy-number output
sv_cn_dtype = [('chr1', '<U50'),
            ('pos1', int),
            ('dir1', '<U1'),
            ('total_copynumber1', int),
            ('no_chrs_bearing_mutation1', int),
            ('chr2', '<U50'),
            ('pos2', int),
            ('dir2', '<U1'),
            ('total_copynumber2', int),
            ('no_chrs_bearing_mutation2', int)]

# SV most likely copy-numbers output
sv_mlcn_dtype =[('chr1', '<U50'),
            ('pos1', int),
            ('dir1', '<U1'),
            ('chr2', '<U50'),
            ('pos2', int),
            ('dir2', '<U1'),
            ('pos1_bb_CN', '<U50'),
            ('pos2_bb_CN', '<U50'),
            ('most_likely_ref_copynumber', int),
            ('most_likely_variant_copynumber', int),
            ('prop_chrs_bearing_mutation', float),
            ('support', int),
            ('depth', int),      
            ('pv', float)]

# SV multiplicities output
sv_mult_dtype =[('chr1', '<U50'),
            ('pos1', int),
            ('dir1', '<U1'),
            ('chr2', '<U50'),
            ('pos2', int),
            ('dir2', '<U1'),
            ('tumour_copynumber', int),
            ('multiplicity', int),
            ('tumour_copynumber_options', '<U50'),
            ('multiplicity_options', '<U50'),
            ('probabilities', '<U50')]

# SNV copy-number output
snv_cn_dtype =  [('chr', '<U50'),
            ('pos', int),
            ('total_copynumber', int),
            ('no_chrs_bearing_mutation', int)]

# SNV most likely copy-number output
snv_mlcn_dtype =[('chr', '<U50'),
            ('pos', int),
            ('bb_CN', '<U50'),
            ('most_likely_ref_copynumber', int),
            ('most_likely_variant_copynumber', int),
            ('prop_chrs_bearing_mutation', float)]

# SNV multiplicities output
snv_mult_dtype =[('chr', '<U50'),
            ('pos', int),
            ('tumour_copynumber', int),
            ('multiplicity', int),
            ('tumour_copynumber_options', '<U50'),
            ('multiplicity_options', '<U50'),
            ('probabilities', '<U50')]
