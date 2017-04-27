import numpy as np

#####################################################################
# Output dtypes
#####################################################################

# SV copy-number output
sv_cn_dtype = [('chr1', 'S50'),
                ('pos1', int),
                ('dir1', 'S1'),
                ('total_copynumber1', int),
                ('no_chrs_bearing_mutation1', int),
                ('chr2', 'S50'),
                ('pos2', int),
                ('dir2', 'S1'),
                ('total_copynumber2', int),
                ('no_chrs_bearing_mutation2', int)]

# SV most likely copy-numbers output
sv_mlcn_dtype =[('chr1', 'S50'),
                ('pos1', int),
                ('dir1', 'S1'),
                ('chr2', 'S50'),
                ('pos2', int),
                ('dir2', 'S1'),
                ('pos1_bb_CN', 'S50'),
                ('pos2_bb_CN', 'S50'),
                ('most_likely_ref_copynumber', int),
                ('most_likely_variant_copynumber', int),
                ('prop_chrs_bearing_mutation', float),
                ('support', int),
                ('depth', int),
                ('bb_var_copynumber_frac', float),
                ('pv', float),
                ('pv_deviance_from_vaf', float)]

# SV multiplicities output
sv_mult_dtype =[('chr1', 'S50'),
                ('pos1', int),
                ('dir1', 'S1'),
                ('chr2', 'S50'),
                ('pos2', int),
                ('dir2', 'S1'),
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
                ('prop_chrs_bearing_mutation', float),
                ('bb_var_copynumber_frac', float),
                ('pv', float),
                ('pv_deviance_from_vaf', float)]

# SNV multiplicities output
snv_mult_dtype =[('chr', 'S50'),
                ('pos', int),
                ('tumour_copynumber', int),
                ('multiplicity', int),
                ('tumour_copynumber_options', 'S50'),
                ('multiplicity_options', 'S50'),
                ('probabilities', 'S50')]
