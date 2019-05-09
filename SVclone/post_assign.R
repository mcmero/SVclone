library(ccube)
library(dplyr)
source('SVclone/write_output.R')

args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0 | args[1] == "-h") {

    cat("
        Perform post-assignment of SVs.

        Usage:
        Rscript post_assign.R <SV RData results> <SNV RData results> <results folder> <sample> [--joint]\n\n

        The --joint flag indicates that variants will be post-assigned to
        a joint SV + SNV model. By default, SVs will be assigned to the
        provided SNV model.
    ")
    q(save="no")
}

svres <- args[1]
snvres <- args[2]
resultFolder <- args[3]
sample <- args[4]
joint_model  <- length(grep("--joint", args)) > 0

if (!file.exists(svres) | !file.exists(snvres)) {
    print('SV or SNV data file do not exist. Exiting.')
    q(save="no")
}

load(svres)
load(snvres)

system(paste('mkdir -p', resultFolder))

svdata <- doubleBreakPtsRes$ssm
snvdata <- snvRes$ssm
if (joint_model) {
    print('Post-assigning using a joint SV + SNV model...')
    postAssignResSVs <- RunPostAssignPipeline(snvRes = snvRes$res,
                                              svRes = doubleBreakPtsRes$res,
                                              mydata = svdata)
    save(postAssignResSVs, file=paste0(resultFolder, sample, "_ccube_postAssign_sv_results.RData"))
    write_sv_output(postAssignResSVs, resultFolder, sample)
    MakeCcubeStdPlot_sv(res = postAssignResSVs$res, ssm = postAssignResSVs$ssm,
                        printPlot = T, fn = paste0(resultFolder, sample, "_ccube_sv_postAssign_results.pdf"))

    system(paste0('mkdir -p ', resultFolder, '/snvs'))
    postAssignResSNVs <- RunPostAssignPipeline(snvRes = snvRes$res,
                                               svRes = doubleBreakPtsRes$res,
                                               mydata = snvdata)
    save(postAssignResSNVs, file=paste0(resultFolder, '/snvs/', sample, "_ccube_postAssign_snv_results.RData"))
    write_snv_output(postAssignResSNVs, paste0(resultFolder, '/snvs/'), sample)
    MakeCcubeStdPlot(res = postAssignResSNVs$res, ssm = postAssignResSNVs$ssm,
                        printPlot = T, fn = paste0(resultFolder, '/snvs/', sample, "_ccube_snv_postAssign_results.pdf"))
} else {
    print('Post-assigning SVs using SNV results...')
    postAssignRes <- RunPostAssignPipeline(svRes = snvRes$res,
                                           mydata = svdata)
    save(postAssignRes, file=paste0(resultFolder, sample, "_ccube_sv_postAssign_results.RData"))
    write_sv_output(postAssignRes, resultFolder, sample)
    MakeCcubeStdPlot_sv(res = postAssignRes$res, ssm = postAssignRes$ssm,
                        printPlot = T, fn = paste0(resultFolder, sample, "_ccube_sv_postAssign_results.pdf"))
}

