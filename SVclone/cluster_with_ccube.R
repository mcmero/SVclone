library(dplyr)
library(ccube)
library(doParallel)
source('SVclone/write_output.R')
options(error = quote({dump.frames(to.file=TRUE); q()}))

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3) {
    args <- c("--help")
}

if("--help" %in% args) {
    cat("
        Cluster provided input file using ccube

        Usage:
        Rscript cluster_with_ccube.R <input_file> <out_dir> <sample_name> --cores=<N> --clusmax=<N> --repeat=<N> --maxiter=<N>\n\n")
    q(save="no")
}

ccfile <- args[1]
outdir <- args[2]
sample <- args[3]
clusmax  <- grep("--clusmax=", args, value=TRUE)
cores  <- grep("--cores=", args, value=TRUE)
numOfRepeat <- grep("--repeat", args, value=TRUE)
maxiter <- grep("--maxiter", args, value=TRUE)

if (length(cores) == 0) {
    cores <- 1
} else {
    cores <- as.numeric(strsplit(cores, "=")[[1]][2])
}
registerDoParallel(cores=cores)

if (length(clusmax) == 0) {
    clusmax <- 6
} else {
    clusmax <- as.numeric(strsplit(clusmax, "=")[[1]][2])
}

if (length(numOfRepeat) == 0) {
    numOfRepeat <- 5
} else {
    numOfRepeat <- as.numeric(strsplit(numOfRepeat, "=")[[1]][2])
}

if (length(maxiter) == 0) {
    maxiter <- 1000
} else {
    maxiter <- as.numeric(strsplit(maxiter, "=")[[1]][2])
}

if (!file.exists(ccfile)) {
    cat("ccube input file does not exist!\n")
    q(save="no")
}

cc_input <- read.delim(ccfile, sep="\t")
if (nrow(cc_input) == 0) {
    cat("ccube input is empty!\n")
    q(save="no")
}

#ensure clusters <= number of data points
clusmax <- min(clusmax, nrow(cc_input))
numOfClusterPool = 1:clusmax

multiCore <- FALSE
if (cores > 1) {multiCore <- TRUE}

#################################################################################################################
# Cluster
#################################################################################################################
resultFolder <- paste(outdir, "ccube_out/", sep="/")
system(paste("mkdir -p", resultFolder))

# clean up any possible NA copynumbers
cn_cols <- grep("cn",colnames(cc_input))
cc_input[,cn_cols][is.na(cc_input[,cn_cols])] <- 0

is_sv_data = "var_counts1" %in% colnames(cc_input)
if (is_sv_data) {
    doubleBreakPtsRes <- RunCcubePipeline(dataFolder = outdir, sampleName = sample,
                                          ssm = cc_input, modelSV = T,
                                          numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                          runAnalysis = T, runQC = T, maxiter = maxiter,
                                          ccubeResultRDataFile = paste(resultFolder, "ccube_sv_results.RData", sep="/"),
                                          multiCore = multiCore, basicFormats = F, allFormats = F, returnAll = T)
    save(doubleBreakPtsRes, file=paste0(resultFolder, sample, "_ccube_sv_results.RData"))
    MakeCcubeStdPlot_sv(res = doubleBreakPtsRes$res, ssm = doubleBreakPtsRes$ssm,
                        printPlot = T, fn = paste0(resultFolder, sample, "_ccube_sv_results.pdf"))
} else {
    resultFolder <- paste(resultFolder, "snvs/", sep="/")
    system(paste("mkdir -p", resultFolder))
    snvRes <- RunCcubePipeline(dataFolder = outdir, sampleName = sample,
                                ssm = cc_input, resultFolder = resultFolder,
                                ccubeResultRDataFile = paste(resultFolder, "ccube_snv_results.RData", sep="/"),
                                numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                runAnalysis = T, runQC = T, multiCore = multiCore,
                                basicFormats = F, allFormats = F, returnAll = T)
    save(snvRes, file=paste0(resultFolder, sample, "_ccube_snv_results.RData"))
    MakeCcubeStdPlot(res = snvRes$res, ssm = snvRes$ssm, printPlot = T, fn = paste0(resultFolder, sample, "_ccube_snv_results.pdf"))
}

#################################################################################################################
# Write output
#################################################################################################################

if (is_sv_data) {
    write_sv_output(doubleBreakPtsRes, resultFolder, sample)
} else {
    write_snv_output(snvRes, resultFolder, sample)
}
