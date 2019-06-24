write_sv_output <- function(doubleBreakPtsRes, resultFolder, sample) {
    res <- doubleBreakPtsRes$res
    uniqLabels <- unique(res$label)
    cc_input <- doubleBreakPtsRes$ssm

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
