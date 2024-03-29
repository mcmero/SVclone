library(dplyr)
library(ccube)
library(doParallel)

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
                                          multiCore = multiCore, basicFormats = F, allFormats = F)
    save(doubleBreakPtsRes, file=paste0(resultFolder, sample, "_ccube_sv_results.RData"))
    if(length(doubleBreakPtsRes$res) > 0) {
        MakeCcubeStdPlot_sv(res = doubleBreakPtsRes$res, ssm = doubleBreakPtsRes$ssm,
                            printPlot = T, fn = paste0(resultFolder, sample, "_ccube_sv_results.pdf"))
    }
} else {
    resultFolder <- paste(resultFolder, "snvs/", sep="/")
    system(paste("mkdir -p", resultFolder))
    snvRes <- RunCcubePipeline(dataFolder = outdir, sampleName = sample,
                                ssm = cc_input, resultFolder = resultFolder,
                                ccubeResultRDataFile = paste(resultFolder, "ccube_snv_results.RData", sep="/"),
                                numOfClusterPool = numOfClusterPool, numOfRepeat = numOfRepeat,
                                runAnalysis = T, runQC = T, multiCore = multiCore,
                                basicFormats = F, allFormats = F)
    save(snvRes, file=paste0(resultFolder, sample, "_ccube_snv_results.RData"))
    MakeCcubeStdPlot(res = snvRes$res, ssm = snvRes$ssm, printPlot = T, fn = paste0(resultFolder, sample, "_ccube_snv_results.pdf"))
}

#################################################################################################################
# Write output
#################################################################################################################

write_sv_output <- function(doubleBreakPtsRes, resultFolder, sample) {
    if(length(doubleBreakPtsRes$res) > 0) {
        res <- doubleBreakPtsRes$res
        uniqLabels <- unique(res$label)
        cc_input <- doubleBreakPtsRes$ssm
    } else if (nrow(doubleBreakPtsRes) == 1) {
        # handle case if clustering has been run with a single SV
        res <- NULL
        res$labels <- 1
        res$full.model$responsibility <- 1
        uniqLabels <- 1
        cc_input <- dplyr::mutate(dplyr::rowwise(doubleBreakPtsRes),
                        ccube_ccf1 = MapVaf2CcfPyClone(vaf1,
                                                       purity,
                                                       normal_cn,
                                                       total_cn1,
                                                       total_cn1,
                                                       1, constraint = F),
                        ccube_ccf2 = MapVaf2CcfPyClone(vaf2,
                                                       purity,
                                                       normal_cn,
                                                       total_cn2,
                                                       total_cn2,
                                                       1, constraint = F))
        res$full.model$ccfMean <- mean(cc_input$ccube_ccf1, cc_input$ccube_ccf2)
    } else {
        stop("Invalid SV output.")
    }

    svs <- do.call(rbind, strsplit(as.character(cc_input$mutation_id), "_", fixed = T))
    sv1 <- do.call(rbind, strsplit(as.character(svs[,1]), ":", fixed = T))
    sv2 <- do.call(rbind, strsplit(as.character(svs[,2]), ":", fixed = T))
    svs <- data.frame(sv1, sv2)
    colnames(svs) <- c("chr1", "pos1", "dir1", "chr2", "pos2", "dir2")

    #### assignment probability ####
    mutAssign <- svs
    if (length(uniqLabels) == 1) {
      mutR = data.frame(res$full.model$responsibility)
      colnames(mutR) <- "cluster0"
    } else {
      mutR <- data.frame(res$full.model$responsibility[, sort(uniqLabels)])
      colnames(mutR) <- paste0("cluster", seq_along(uniqLabels)-1)
    }
    mutAssign <- data.frame(mutAssign, mutR)
    fn <- paste0(resultFolder, "/", sample, "_assignment_probability_table.txt")
    write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)

    #### multiplicity ####
    mult <- svs
    mult$tumour_copynumber1 <- cc_input$total_cn1
    mult$tumour_copynumber2 <- cc_input$total_cn2
    mult$multiplicity1 <- cc_input$ccube_mult1
    mult$multiplicity2 <- cc_input$ccube_mult2

    fn <- paste0(resultFolder, "/", sample, "_multiplicity.txt")
    write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)

    #### subclonal structure ####
    cellularity <- unique(cc_input$purity)
    subc <- as.data.frame(table(res$label), stringsAsFactors = F)
    subc <- dplyr::rename(subc, cluster = Var1, n_ssms = Freq)
    subc$proportion <- res$full.model$ccfMean[as.integer(subc$cluster)] * cellularity
    fn <- paste0(resultFolder, "/", sample, "_subclonal_structure.txt")
    write.table(subc, file = fn, sep = "\t", row.names = F, quote = F)

    #### cluster certainty ####
    ccert <- svs
    ccert$most_likely_assignment <- res$label
    ccert$average_proportion1 <- cc_input$ccube_ccf1 * cellularity
    ccert$average_proportion2 <- cc_input$ccube_ccf2 * cellularity
    fn <- paste0(resultFolder, "/", sample, "_cluster_certainty.txt")
    write.table(ccert, file = fn, sep = "\t", row.names = F, quote = F)
}

write_snv_output <- function(snvRes, resultFolder, sample) {
    res <- snvRes$res
    uniqLabels <- unique(res$label)
    cc_input <- snvRes$ssm

    snvs <- data.frame(do.call(rbind, strsplit(as.character(cc_input$mutation_id), "_", fixed = T)))
    colnames(snvs) <- c("chr", "pos")

    #### assignment probability ####
    mutAssign <- snvs
    if (length(uniqLabels) == 1) {
      mutR = data.frame(res$full.model$responsibility)
      colnames(mutR) <- "cluster0"
    } else {
      mutR <- data.frame(res$full.model$responsibility[, sort(uniqLabels)])
      colnames(mutR) <- paste0("cluster", seq_along(uniqLabels)-1)
    }
    mutAssign <- data.frame(mutAssign, mutR)
    fn <- paste0(resultFolder, "/", sample, "_assignment_probability_table.txt")
    write.table(mutAssign, file = fn, sep = "\t", row.names = F, quote = F)

    #### multiplicity ####
    mult <- snvs
    mult$tumour_copynumber <- cc_input$total_cn
    mult$multiplicity <- cc_input$ccube_mult

    fn <- paste0(resultFolder, "/", sample, "_multiplicity.txt")
    write.table(mult, file = fn, sep = "\t", row.names = F, quote = F)

    #### subclonal structure ####
    cellularity <- unique(cc_input$purity)
    subc <- as.data.frame(table(res$label), stringsAsFactors = F)
    subc <- dplyr::rename(subc, cluster = Var1, n_ssms = Freq)
    subc$proportion <- res$full.model$ccfMean[as.integer(subc$cluster)] * cellularity
    fn <- paste0(resultFolder, "/", sample, "_subclonal_structure.txt")
    write.table(subc, file = fn, sep = "\t", row.names = F, quote = F)

    #### cluster certainty ####
    ccert <- snvs
    ccert$most_likely_assignment <- res$label
    ccert$average_proportion <- cc_input$ccube_ccf * cellularity
    fn <- paste0(resultFolder, "/", sample, "_cluster_certainty.txt")
    write.table(ccert, file = fn, sep = "\t", row.names = F, quote = F)
}

if (is_sv_data) {
    write_sv_output(doubleBreakPtsRes, resultFolder, sample)
} else {
    write_snv_output(snvRes, resultFolder, sample)
}
