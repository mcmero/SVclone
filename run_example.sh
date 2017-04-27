#!/bin/sh

input=example_data/tumour_p80_DEL_svs_simple.txt
bam=example_data/tumour_p80_DEL_sv_extract_sorted.bam
sample=tumour_p80_DEL

./SVclone.py annotate -i $input -b $bam -s $sample --sv_format simple -cfg svclone_test.ini

./SVclone.py count -i ${sample}/${sample}_svin.txt -b $bam -s $sample -cfg svclone_test.ini

./SVclone.py filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt -cfg svclone_test.ini

./SVclone.py cluster -s $sample -cfg svclone_test.ini

Rscript post_process_fit_diagnostics.R $sample $sample

# uncomment below to test run with SNVs:

#./SVclone.py filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt --snvs example_data/tumour_p80_DEL_snvs.vcf --snv_format mutect

#./SVclone.py cluster -s $sample

#Rscript post_process_fit_diagnostics.R $sample $sample --snvs
