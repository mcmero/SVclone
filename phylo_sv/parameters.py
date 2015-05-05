import numpy as np

# constant parameters

tr                      = 5    # threshold by how much read has to overlap breakpoint
window                  = 500  # base-pairs considered to the left and right of the break
subclone_threshold      = 0.05 # throw out any subclones with frequency lower than this value
subclone_diff           = 0.10 # merge any clusters within this range

# parameters extracted for each read from BAMs
read_dtype = [('query_name', 'S150'), ('chrom', 'S50'), ('ref_start', int), ('ref_end', int), \
              ('align_start', int), ('align_end', int), ('len', int), ('ins_len', int), ('is_reverse', np.bool)]

