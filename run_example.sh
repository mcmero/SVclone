#!/bin/sh

sample=tumour_p80_DEL
input=example_data/${sample}_svs_simple.txt
bam=example_data/${sample}_sv_extract_sorted.bam

echo 'Annotating breakpoints...'
svclone annotate -i $input -b $bam -s $sample --sv_format simple -cfg svclone_test.ini

echo 'Counting breakpoint reads...'
svclone count -i ${sample}/${sample}_svin.txt -b $bam -s $sample

echo 'Filtering out low-confidence breaks...'
svclone filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt

echo 'Clustering...'
svclone cluster -s $sample

if [ -f tumour_p80_DEL/ccube_out/tumour_p80_DEL_ccube_sv_results.RData ] ; then
    echo 'Successfully ran test sample to completion!'
fi
