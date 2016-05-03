#!/bin/sh

input=example_data/tumour_p80_DEL_svs_simple.txt
bam=example_data/tumour_p80_DEL_sv_extract_sorted.bam
sample=tumour_p80_DEL

./SVclone.py identify -i $input -b $bam -s $sample --sv_format simple

./SVclone.py count -i ${sample}/${sample}_svin.txt -b $bam -s $sample

./SVclone.py filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt

./SVclone.py cluster -s $sample

Rscript post_process_fit_diagnostics.R $sample $sample
