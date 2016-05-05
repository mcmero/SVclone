import numpy as np

#####################################################################
# Data dtypes
#####################################################################

bp_dtype = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S1')]

sv_dtype = [('ID', 'int64'),
            ('chr1', 'S20'),
            ('pos1', 'int64'),
            ('dir1', 'S1'),
            ('chr2', 'S20'),
            ('pos2', 'int64'),
            ('dir2', 'S1'),
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
            ('chr1', 'S20'),
            ('pos1', 'int64'),
            ('dir1', 'S1'),
            ('chr2', 'S20'),
            ('pos2', 'int64'),
            ('dir2', 'S1'),
            ('classification', 'S100'),
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
