#!/bin/sh

input=example_data/tumour_p80_DEL_svs_simple.txt
bam=example_data/tumour_p80_DEL_sv_extract_sorted.bam
sample=tumour_p80_DEL

echo 'Annotating breakpoints...'
python SVclone.py annotate -i $input -b $bam -s $sample --sv_format simple -cfg svclone_test.ini

echo 'Counting breakpoint reads...'
python SVclone.py count -i ${sample}/${sample}_svin.txt -b $bam -s $sample -cfg svclone_test.ini

echo 'Filtering out low-confidence breaks...'
python SVclone.py filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt -cfg svclone_test.ini

python SVclone.py cluster -s $sample -cfg svclone_test.ini

echo 'Post-assigning variants...'
python SVclone.py post_assign -s $sample --svs ${sample}/${sample}_svinfo.txt -cfg svclone_test.ini

echo 'Drawing diagnostic plots...'
Rscript post_process_fit_diagnostics.R $sample $sample --map

if [[ -f tumour_p80_DEL/best_run_svs_post_assign/tumour_p80_DEL_best_run_svs_post_assign_best_fit.pdf ]] 
    then
        echo 'Successfully ran test sample to completion!'
fi

# uncomment below to test run with SNVs:
#echo 'Filtering out low-confidence breaks, adding SNVs...'
#python SVclone.py filter -s $sample -i ${sample}/${sample}_svinfo.txt -p example_data/purity_ploidy.txt --snvs example_data/tumour_p80_DEL_snvs.vcf --snv_format mutect

#python SVclone.py cluster -s $sample -cfg svclone_test.ini

#echo 'Post-assigning variants...'
#python SVclone.py post_assign -s $sample --svs ${sample}/${sample}_svinfo.txt -cfg svclone_test.ini --snvs example_data/tumour_p80_DEL_    snvs.vcf --snv_format mutect

#echo 'Drawing diagnostic plots...'
#Rscript post_process_fit_diagnostics.R $sample $sample --colcus --map
