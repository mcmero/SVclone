args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 3 | args[1]=='-h') {
    print('Usage: Rscript <workingdir> <sample id> <run> --snvs. To select best run, use best as run')
    quit()
}

library('mcclust')
setwd(args[1])
id <- args[2]
run <- args[3]
snvs <- FALSE

if(length(args) > 2) {
    snvs <- '--snvs' %in% args
}

if (run == 'best') {
    if (snvs) run <- '/best_run_snvs'
    else run <- '/best_run_svs'
}

z_trace <- paste(id, run, '/z_trace.txt.gz', sep='')
z <- read.table(gzfile(z_trace), header=F)

psm <- comp.psm(as.matrix(z)+1)

outfile <- paste(id, run, '/coclustering_matrix.txt', sep='')
write.table(psm, file = outfile, quote = F, sep = '\t', col.names = F, row.names = F)
