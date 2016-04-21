#!/bin/sh

input=example_data/tumour_p80_DEL_svs_simple.txt
bam=example_data/tumour_p80_DEL_sv_extract_sorted.bam
sample=tumour_p80_DEL

./SVClone.py identify -i $input -b $bam -s $sample --simple

./SVClone.py count -i ${sample}/${sample}_svin.txt -b $bam -s $sample

./SVClone.py filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt

./SVClone.py cluster -s $sample -n 4 --map --burn 2000

Rscript post_process_fit_diagnostics.R $sample $sample --map
