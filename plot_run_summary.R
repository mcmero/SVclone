#plot run summary

args <- commandArgs(trailingOnly = TRUE)
if(length(args)==0 | args[1]=='-h') {
    print('Usage: Rscript <wd> <p1> <p2> <postfix> --snvs')
}

source('~/Dropbox/bioinfo/svclone/001_analysis/mix_output_functions.R')
wd <- paste(args[1], '/', sep='')
p1 <- as.numeric(args[2])
p2 <- as.numeric(args[3])
ps <- args[4]
snvs <- FALSE
if(is.na(ps) | ps=='--snvs'){ps <- ''}
if('--snvs' %in% args){snvs <- TRUE}

if (ps != '') {ps <- paste('_', ps, sep='')}

wd_truth <- '~/Dropbox/bioinfo/svclone/001/'
if (snvs) {
    wd_truth <- "~/Dropbox/bioinfo/svclone/001_snvs/snv_counts_mutec_ss/"
}

mixes <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)
mixes <- data.frame(bM=mixes, gM=100 - mixes)

p_name <- paste('m', p1, p2, sep='')
c_name <- paste('i', p1, p2, sep='')
p <- plot_mix(wd, wd_truth, p1, p2, ps, snvs)
pi <- plot_ics(wd, p1, p2, ps, snvs)

base_name <- paste('001bM_p', p1/100, '_001gM_p', p2/100, sep='')
pdf_name <- paste(wd, '/', base_name, ps, '/', base_name, '_runs.pdf', sep='')
pdf(pdf_name, height = 6, width = 16)
grid.arrange(p, pi, ncol = 2)
dev.off()
