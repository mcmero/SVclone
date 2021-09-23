#!/bin/sh

sample=tumour_p80_DEL
input=example_data/${sample}_svs_simple.txt
bam=example_data/${sample}_sv_extract_sorted.bam

echo 'Annotating breakpoints...'
svclone annotate -i $input -b $bam -s $sample --sv_format simple -cfg svclone_test.ini

echo 'Counting breakpoint reads...'
svclone count -i ${sample}/${sample}_svin.txt -b $bam -s $sample

echo 'Filtering out low-confidence breaks, adding SNVs...'
svclone filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt --snvs example_data/tumour_p80_DEL_snvs.vcf --snv_format mutect

echo 'Clustering...'
svclone cluster -s $sample

echo 'Post-assigning variants...'
svclone postassign -s $sample --joint

if [ -f ${sample}/ccube_out/post_assign/snvs/${sample}_ccube_postAssign_snv_results.RData ]
    then
        echo 'Successfully ran test sample to completion!'
fi
