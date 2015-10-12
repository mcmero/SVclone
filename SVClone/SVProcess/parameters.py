import numpy as np

# PREPROCESSING PARAMETERS
tr              = 6    # "wobble length" tolerance threshold which we allow breaks to be inexact
window          = 500  # base-pairs considered to the left and right of the break
norm_overlap    = 10   # minimum basepairs a "normal" read must overlap break to be counted
bp_dtype        = [('chrom','S20'),('start', int), ('end', int), ('dir', 'S1')]

# parameters extracted for each read from BAMs
read_dtype      = [('query_name',       'S150'), 
                   ('chrom',            'S50'), 
                   ('ref_start',        'int64'), 
                   ('ref_end',          'int64'), 
                   ('align_start',      'int64'), 
                   ('align_end',        'int64'), 
                   ('len',              'int64'), 
                   ('ins_len',          'int64'), 
                   ('is_reverse',       np.bool)]

sv_dtype        = [('ID',            'int64'),
                   ('bp1_chr',          'S20'),
                   ('bp1_pos',          'int64'),
                   ('bp1_dir',          'S1'),
                   ('bp2_chr',          'S20'),
                   ('bp2_pos',          'int64'),
                   ('bp2_dir',          'S1'),
                   ('classification',   'S100')] 

# dtypes for SV output file
sv_out_dtype    = [('ID',               'int64'),
                   ('bp1_chr',          'S20'),
                   ('bp1_pos',          'int64'),
                   ('bp1_dir',          'S1'),
                   ('bp2_chr',          'S20'),
                   ('bp2_pos',          'int64'),
                   ('bp2_dir',          'S1'),
                   ('classification',   'S100'),
                   ('bp1_split_norm',   'int64'),
                   ('bp1_norm_olap_bp', 'int64'),
                   ('bp1_span_norm',    'int64'),
                   ('bp1_win_norm',     'int64'), 
                   ('bp1_split',        'int64'), 
                   ('bp1_sc_bases',     'int64'), 
                   ('bp1_total_reads',  'int64'),
                   ('bp2_split_norm',   'int64'), 
                   ('bp2_norm_olap_bp', 'int64'),
                   ('bp2_span_norm',    'int64'), 
                   ('bp2_win_norm',     'int64'), 
                   ('bp2_split',        'int64'), 
                   ('bp2_sc_bases',     'int64'), 
                   ('bp2_total_reads',  'int64'),
                   ('spanning',         'int64'),
                   ('norm1',            'int64'),
                   ('norm2',            'int64'), 
                   ('support',          'int64'),
                   ('vaf1',             float),
                   ('vaf2',             float)]

# Socrates fields
bp1_pos         = 'C1_anchor'
bp1_dir         = 'C1_anchor_dir'
bp2_pos         = 'C2_anchor'
bp2_dir         = 'C2_anchor_dir'
avg_mapq1       = 'C1_avg_realign_mapq'
avg_mapq2       = 'C2_avg_realign_mapq'
repeat1         = 'repeat1'
repeat2         = 'repeat2'
