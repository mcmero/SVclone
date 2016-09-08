# get a metric for stability of clusters

args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0 | args[1]=='-h') {
    print('Usage: Rscript <workingdir> <sample id> --map --clus_stability --snvs --sc')
}

setwd(args[1])
id <- args[2]
map <- FALSE
clus_stab <- FALSE
snvs <- FALSE
handle_sc <- FALSE
allowed_ics <- c( 'AICc', 'AIC', 'BIC', 'svc_IC')

if(length(args) > 2) {
    map <- '--map' %in% args
    clus_stab <- '--clus_stability' %in% args
    snvs <- '--snvs' %in% args
    handle_sc <- '--sc' %in% args
}

snv_dir <- ''
if (snvs) {
    snv_dir <- 'snvs/'
}

require(combinat)
require(RColorBrewer)
require(gplots)
require(ggplot2)
require(gridExtra)
require(reshape)
require(gtools)

# Get number of runs
runs <- c()
for (dir in list.dirs()) {
    cur_dir <- strsplit(dir,'/')[[1]]
    if (length(cur_dir) <= 1)
        next
    cur_dir <- cur_dir[2]
    if (substring(cur_dir,1,3) == 'run') {
        runs <- c(runs,cur_dir)
    }
}
runs <- unique(runs)

############################################################################################
# Cluster stability of runs
############################################################################################

compare_runs <- function(run1, run2) {
    clus_cert_file <- paste(run1, '/', snv_dir, id, '_cluster_certainty.txt', sep='')
    clus_cert1 <- read.delim(clus_cert_file, sep='\t', stringsAsFactors=F)

    clus_cert_file <- paste(run2, '/', snv_dir, id, '_cluster_certainty.txt', sep='')
    clus_cert2 <- read.delim(clus_cert_file, sep='\t', stringsAsFactors=F)

    clus1 <- as.numeric(names(table(clus_cert1$most_likely_assignment)))
    dist <- 0
    for (clus in clus1) {
        if (sum(clus_cert1$most_likely_assignment == clus) <= 1)
            next

        cert1_tmp <- clus_cert1[clus_cert1$most_likely_assignment == clus,]
        cert2_tmp <- clus_cert2[clus_cert1$most_likely_assignment == clus,]

        points_in_clus <- nrow(cert1_tmp)
        combs <- t(combn(points_in_clus, 2))
        for (i in 1:nrow(combs)) {
            pair <- combs[i,]
            clus_in_2 <- cert2_tmp[pair, 'most_likely_assignment']
            if (clus_in_2[1] != clus_in_2[2]) {
                dist <- dist + 1
            }
        }
    }

    total_combs <- ncol(combn(nrow(clus_cert1), 2))
    #print(paste('Dist:', dist, ';', 'combs:', total_combs))
    metric <- 1 - (dist / total_combs)
    #print('-------------------')
    return(metric)
}

if (clus_stab) {
    print('Comparing clustering stability between runs...')
    all_comps <- combn(runs,2)
    n <- length(runs)

    results <- matrix(0, n, n)
    colnames(results) <- runs
    rownames(results) <- runs

    for (i in 1:nrow(results)) {
        for (j in 1:ncol(results)) {
            run1 <- rownames(results)[i]
            run2 <- colnames(results)[j]
            if (run1 == run2) {
                results[i, j] <- 1.
            } else {
                results[i, j] <- compare_runs(run1, run2)
            }
        }
    }

    # cluster stability metric plot
    pdf(paste(id, '_cluster_stability_heatmap.pdf',sep=''),height=6)
    cols <- colorRampPalette(brewer.pal(9,'Blues'))(30)
    heatmap.2(results, trace='none', Rowv=F, Colv=F, dendrogram='none', col=cols)
    dev.off()

    results <- cbind(rownames(results), results)
    write.table(results, paste(id,'_cluster_stability.csv', sep=''), sep=',', quote=F, row.names=F)
}

# cluster proportions plot
clusts <- NULL
for (run in runs) {
    sv_clust <- read.table(paste(run, '/', snv_dir, id, '_subclonal_structure.txt', sep=''),
                           header=T, sep='\t', stringsAsFactors=F)
    clusts <- rbind(clusts, data.frame(run=run, n_ssms=sv_clust$n_ssms, cluster=sv_clust$cluster))
}

pdf(paste(id, '_cluster_hist.pdf',sep=''),height=4, width=max(3,length(runs)*0.6))
ggplot(clusts, aes(x=factor(run), y=n_ssms, fill=factor(cluster))) + geom_bar(stat='identity')
dev.off()

############################################################################################
# Plot AIC and BIC for runs
############################################################################################

if (length(args)>2 & map) {
    print('Plotting AIC & BIC metrics...')
    ic_table <- NULL
    for (run in runs) {
        ic <- read.table(paste(run, '/', snv_dir, id, '_fit.txt', sep=''), sep='\t', header=F, stringsAsFactors = F)
        ic <- ic[ic$V1%in%allowed_ics,]
        ic <- cbind(run=run, ic)
        ic_table <- rbind(ic_table, ic)
    }

    ic_tmp <- ic_table[ic_table$V1!='svc_IC',]
    ic_tmp$run <- factor(ic_tmp$run, levels=mixedsort(unique(as.character(ic_tmp$run))))
    ic_plot <- ggplot(ic_tmp, aes(y=V2, x=run, group=V1, color=factor(V1))) + ylab('value') + geom_line()
    pdf(paste(id, 'aic_bic_plot.pdf',sep='_'), height=4)
    print(ic_plot)
    dev.off()

    # plot SVclone's metric
    ic_tmp <- ic_table[ic_table$V1=='svc_IC',]
    ic_tmp$run <- factor(ic_tmp$run, levels=mixedsort(unique(as.character(ic_tmp$run))))
    svic_plot <- ggplot(ic_tmp, aes(y=V2, x=run, group=V1, color=factor(V1))) + ylab('value') + geom_line()
    pdf(paste(id, 'svc_ic_plot.pdf',sep='_'), height=4)
    print(svic_plot)
    dev.off()

    ic_table <- cast(ic_table, run~V1, value='V2')
    min_bic <- ic_table[min(ic_table$BIC)==ic_table$BIC,]
    min_bic$AIC <- min_bic$run
    min_bic$run <- 'min_BIC'
    min_bic$AICc <- min_bic$BIC
    min_bic$BIC <- NA; min_bic$svc_IC <- NA

    min_aic <- ic_table[min(ic_table$AIC)==ic_table$AIC,]
    min_aic$AICc <- min_aic$AIC
    min_aic$AIC <- min_aic$run
    min_aic$run <- 'min_AIC'
    min_aic$BIC <- NA; min_aic$svc_IC <- NA

    min_aicc <- ic_table[min(ic_table$AICc)==ic_table$AICc,]
    min_aicc$AIC <- min_aicc$run
    min_aicc$run <- 'min_AICc'
    min_aicc$BIC <- NA; min_aicc$svc_IC <- NA

    min_svcic <- ic_table[min(ic_table$svc_IC,na.rm=T)==ic_table$svc_IC,]
    min_svcic$AIC <- min_svcic$run
    min_svcic$run <- 'min_svc_IC'
    min_svcic$AICc <- min_svcic$svc_IC
    min_svcic$BIC <- NA; min_svcic$svc_IC <- NA

    ic_table <- rbind(ic_table, min_bic)
    ic_table <- rbind(ic_table, min_aic)
    ic_table <- rbind(ic_table, min_aicc)
    ic_table <- rbind(ic_table, min_svcic)

    write.table(ic_table, paste(id,'_aic_bic_metrics.csv', sep=''), sep=',', quote=F, row.names=F, na = '')
}

############################################################################################
# Plot histogram of clusters + QQ plots
############################################################################################

get_adjust_factor <- function(dat, pur) {
    cn_v <- dat$most_likely_variant_copynumber
    cn_r <- dat$most_likely_ref_copynumber
    mut_prop <- dat$prop_chrs_bearing_mutation
    vaf <- dat$adjusted_vaf
    frac <- dat$cn_frac
    prob <- (cn_v * mut_prop * pur) / (2 * (1 - pur) + cn_v * pur)
    if (handle_sc) {
        prob <- (cn_v * mut_prop * pur) / (2 * (1 - pur) + (cn_v * pur * frac) + (cn_r * pur * (1. - frac)))
    }
    return((1 / prob))
}

var_file <- paste(id, '_filtered_svs.tsv', sep='')
if (snvs) {
    var_file <- paste(id, '_filtered_snvs.tsv', sep='')
}
svs <- read.table(var_file, header=T, sep='\t', stringsAsFactors=F)
pur <- read.table('purity_ploidy.txt',header=T, sep='\t', stringsAsFactors=F)$purity

gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

ggQQ <- function(dat) {
    p <- ggplot(dat) +
        stat_qq(aes(sample=CCF, colour = factor(most_likely_assignment)), alpha = 0.5)

    dat_tmp <- dat[!is.na(dat$CCF),]

    clusts <- as.numeric(names(table(dat_tmp$most_likely_assignment)))
    cols <- gg_color_hue(length(clusts))

    for (i in 1:length(clusts)) {
        clus <- clusts[i]
        tmp <- dat_tmp[dat_tmp$most_likely_assignment == clus, 'CCF']
        y <- quantile(tmp, c(0.25, 0.75))
        x <- qnorm(c(0.25, 0.75))
        slope <- diff(y)/diff(x)
        intercept <- y[1L] - slope * x[1L]

        p <- p + geom_abline(slope = slope, intercept = intercept, color=cols[i], alpha=0.5)
    }

    return(p)
}

get_frac <- function(x, snvs) {
    variant_cn <- x['most_likely_variant_copynumber']
    gtype <- x['gtype']
    if (!snvs) {
        side <- x['preferred_side']
        gtype <- x['gtype1']
        if(side==1) {
            gtype <- x['gtype2']
        }
    }
    sc <- strsplit(gtype, '\\|')[[1]]
    sc <- strsplit(sc,',')
    sc1 <- as.numeric(sc[[1]])
    sc1 <- c(sum(sc1[1:2]),sc1[3])
    if(length(sc) > 1) {
        sc2 <- as.numeric(sc[[2]])
        sc2 <- c(sum(sc2[1:2]),sc2[3])
        if (sc1[1] == variant_cn) {
            return(sc1[2])
        } else {
            return(sc2[2])
        }
    }
    return(sc1[2])
}

if (map) {
    ic_table <- NULL
    for (run in runs) {
        ic <- read.table(paste(run, '/', snv_dir, id, '_fit.txt', sep=''), sep='\t', header=F, stringsAsFactors = F)
        ic <- ic[ic$V1%in%allowed_ics,]
        ic <- cbind(run=run, ic)
        ic_table <- rbind(ic_table, ic)
    }
    ic_table$rank <- NULL
    for (aic in allowed_ics) {
        ic_table[ic_table$V1 == aic, 'rank'] <-
            rank(ic_table[ic_table$V1 == aic, 'V2'])
    }
}

for (run in runs) {
    mlcn <- read.table(paste(run, '/', snv_dir, id, '_most_likely_copynumbers.txt', sep=''), header=T, sep='\t', stringsAsFactors=F)
    dat <- NULL
    if (snvs) {
        dat <- merge(svs, mlcn, by.x=c(1,2), by.y=c(1,2))
        certain <- read.table(paste(run, '/', snv_dir, id, '_cluster_certainty.txt', sep=''), sep='\t', header=T)
        dat <- merge(dat, certain,by.x=c(1,2), by.y=c(1,2))
        dat$adjusted_vaf <- dat$var / (dat$ref + dat$var)
    } else {
        merge_cols <- c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2')
        dat <- merge(svs, mlcn, by.x=merge_cols, by.y=merge_cols)
        certain <- read.table(paste(run, '/', id, '_cluster_certainty.txt', sep=''), sep='\t', header=T)
        dat <- merge(dat, certain,by.x=merge_cols, by.y=merge_cols)
    }

    dat$cn_frac <- apply(dat, 1, function(x){get_frac(x, snvs)})
    dat <- cbind(dat, CCF=get_adjust_factor(dat, pur) * dat$adjusted_vaf)
    dat$CCF <- sapply(dat$CCF,function(x){min(2,x)})

    sv_clust <- read.table(paste(run, '/', snv_dir, id, '_subclonal_structure.txt', sep=''), header=T, sep='\t', stringsAsFactors=F)
#     sv_clust <- sv_clust[sv_clust$n_ssms>1, ]
    dat <- dat[dat$most_likely_assignment%in%sv_clust$cluster,]

    above_ssm_th <- sv_clust$n_ssms / (sum(sv_clust$n_ssms)) >= 0.05
    below_ssm_th <- sv_clust$n_ssms / (sum(sv_clust$n_ssms)) < 0.05
    clus_intercepts <- 1 / pur * as.numeric(sv_clust$proportion[above_ssm_th & sv_clust$n_ssms > 2])
    clus_intercepts_minor <- 1 / pur * as.numeric(sv_clust$proportion[below_ssm_th | sv_clust$n_ssms<=2])

    plot1 <- ggplot(dat, aes(x=as.numeric(dat$CCF),
                    fill=factor(most_likely_assignment), color=factor(most_likely_assignment))) +
                    xlim(0,2) + geom_histogram(alpha=0.3,position='identity',binwidth=0.05)+xlab('CCF') +
                    geom_vline(xintercept=clus_intercepts, colour='blue', size=1)

    if (length(clus_intercepts_minor) > 0) {
        plot1 <- plot1 + geom_vline(xintercept=clus_intercepts_minor,colour='red',lty=2)
    }

    plot2 <- ggplot(dat, aes(x=CCF, y=adjusted_vaf,
                             fill=factor(most_likely_assignment),
                             colour=factor(most_likely_assignment))) +
                    geom_point(size=1) + xlim(0,2) + ylim(0,1) + ylab('VAF')

    # attach table for convenience, also add BIC/AIC
    tabout <- sv_clust[order(as.numeric(sv_clust[,3]),decreasing=TRUE),]
    sc_tab <- tableGrob(tabout, rows=c())
    if (map) {
        ic <- ic_table[ic_table$run==run,]
        rownames(ic) <- ic$V1
        ic <- ic[,-c(1:2)]
        colnames(ic)[1] <- 'value'
        ic$value <- round(ic$value, 4)
        ic_tab <- tableGrob(ic)
        height <- 7+round(nrow(tabout)*0.2)
        pdf(paste(id, run, 'fit.pdf',sep='_'), height=height)
        grid.arrange(arrangeGrob(sc_tab, ic_tab, nrow=1), plot1, plot2, ncol=1)
        dev.off()
    } else {
        height <- 7+round(nrow(tabout)*0.2)
        pdf(paste(id, run, 'fit.pdf',sep='_'), height=height)
        grid.arrange(sc_tab, plot1, plot2, ncol=1)
        dev.off()
    }
}

get_run_info <- function(id, run) {
    scs_file <- paste(run, '/', snv_dir, id, '_subclonal_structure.txt', sep='')
    scs <- read.table(scs_file, sep = '\t', header = T)
    scs <- scs[scs$n_ssms>3,]
    scs <- scs[order(scs$CCF, decreasing=T), ]
    scs$new_cluster <- 1:nrow(scs)

    sv_df <- NULL
    if (snvs) {
        sv_df <- read.table(paste(id, '_filtered_snvs.tsv', sep=''),
                            header=T, sep='\t', stringsAsFactors=F)
    } else {
        sv_df <- read.table(paste(id, '_filtered_svs.tsv', sep=''),
                            header=T, sep='\t', stringsAsFactors=F)
    }
    pur <- read.table(paste('purity_ploidy.txt', sep=''),
                      header=T, sep='\t', stringsAsFactors=F)$purity

    cc_file <-  paste(run, '/', snv_dir, id, '_cluster_certainty.txt', sep = '')
    cc <- read.table(cc_file, sep = '\t', stringsAsFactors = F, header = T)

    mlcn_file <- paste(run, '/', snv_dir, id, '_most_likely_copynumbers.txt', sep='')
    mlcn <- read.table(mlcn_file, header = T, sep = '\t', stringsAsFactors = F)

    merge_cols <- c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2')
    if (snvs) {
        colnames(sv_df)[1] <- 'chr'
        colnames(cc)[1] <- 'chr'
        merge_cols <- c('chr', 'pos')
    }
    dat <- merge(sv_df, mlcn, by.x=merge_cols, by.y=merge_cols)
    dat <- merge(dat, cc, by.x=merge_cols, by.y=merge_cols)
    if (snvs) {
        dat$adjusted_vaf <- dat$var / (dat$ref + dat$var)
    }
    dat$cn_frac <- apply(dat, 1, function(x){get_frac(x, snvs)})
    dat <- cbind(dat, CCF=get_adjust_factor(dat, pur) * dat$adjusted_vaf)

    dat$cluster <- NA
    for (j in 1:nrow(scs)) {
        dat[dat$most_likely_assignment==scs$cluster[j], 'cluster'] <- scs$new_cluster[j]
    }
    if (snvs) {
        dat$sv <- paste(dat$chr, dat$pos, sep='_')
    } else {
        dat$sv <- paste(dat$chr1, dat$pos1, dat$chr2, dat$pos2, sep=':')
    }

    scs$cluster <- scs$new_cluster
    scs <- scs[,1:4]
    scs$variant_proportion <- scs$n_ssms/sum(scs$n_ssms)
    return(list(dat,scs))
}

var_id <- c('chr1','pos1','chr2','pos2')
if (snvs) {var_id <- c('chr', 'pos')}

all_runs_ccfs <- NULL
all_runs_scs <- NULL
for(run in runs) {
    runinf <- get_run_info(id, run)
    dat <- runinf[[1]]
    scs <- runinf[[2]]

    scs$run <- run
    ccfs <- data.frame(sv = apply(dat[,var_id], 1, paste, collapse=':'),
                       CCF = dat$CCF, cluster = dat$cluster,
                       run = run)
    all_runs_ccfs <- rbind(all_runs_ccfs, ccfs)
    all_runs_scs <- rbind(all_runs_scs, scs)
}
all_runs_ccfs$run <- factor(all_runs_ccfs$run, levels=mixedsort(unique(as.character(all_runs_ccfs$run))))

rp <- ggplot(data = all_runs_scs) +
    scale_y_continuous(breaks = seq(0, 2, 0.2), limits = c(0,2)) + xlab('') +
    geom_line(data = all_runs_ccfs, aes(x = run, y = CCF, group = sv,
                                        colour = factor(cluster)), alpha=0.2) +
    geom_jitter(data = all_runs_ccfs, aes(x = run, y = CCF, colour = factor(cluster)),
                position=position_jitter(width=0.4), alpha = 0.5, size = 2) +
    scale_size(range = c(4,25), guide = FALSE) +
    geom_point(aes(x = run, y = CCF, size=variant_proportion), colour = '#4d4d4d', alpha=0.8)

if (map) {
    pdf_name <- paste(id, '_run_summary.pdf', sep='')
    pdf(pdf_name, height = 6, width = 20 * (length(runs)/12))
    grid.arrange(rp, ic_plot, svic_plot, nrow = 1)
    dev.off()
} else {
    pdf_name <- paste(id, '_run_summary.pdf', sep='')
    pdf(pdf_name, height = 6, width = 8 * (length(runs)/10))
    print(rp)
    dev.off()
}
