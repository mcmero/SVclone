import numpy as np

# PREPROCESSING PARAMETERS
tr              = 5    # "wobble length" tolerance threshold which we allow breaks to be inexact
window          = 500  # base-pairs considered to the left and right of the break

# parameters extracted for each read from BAMs
read_dtype      = [('query_name',       'S150'), 
                   ('chrom',            'S50'), 
                   ('ref_start',        int), 
                   ('ref_end',          int), 
                   ('align_start',      int), 
                   ('align_end',        int), 
                   ('len',              int), 
                   ('ins_len',          int), 
                   ('is_reverse',       np.bool)]

# dtypes for SV input file
sv_dtype        = [('bp1_chr',          'S20'),
                   ('bp1_pos',          int),
                   ('bp2_chr',          'S20'), 
                   ('bp2_pos',          int),
                   ('classification',   'S100')]     

# dtypes for SV output file
sv_out_dtype    = [('bp1_split_norm',   int), 
                   ('bp1_span_norm',    int),
                   ('bp1_win_norm',     int), 
                   ('bp1_split',        int), 
                   ('bp1_sc_bases',     int), 
                   ('bp2_split_norm',   int), 
                   ('bp2_span_norm',    int), 
                   ('bp2_win_norm',     int), 
                   ('bp2_split',        int), 
                   ('bp2_sc_bases',     int), 
                   ('spanning',         int),
                   ('norm1',            int),
                   ('norm2',            int), 
                   ('support',          int),
                   ('vaf1',             float),
                   ('vaf2',             float)]
