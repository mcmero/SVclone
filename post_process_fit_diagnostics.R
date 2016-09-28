# get a metric for stability of clusters

args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0 | args[1]=='-h') {
    print('Usage: Rscript <workingdir> <sample id> --map --snvs')
}

setwd(args[1])
id <- args[2]
map <- FALSE
clus_stab <- FALSE
snvs <- FALSE
coclus <- FALSE
allowed_ics <- c( 'AICc', 'AIC', 'BIC', 'svc_IC')

if(length(args) > 2) {
    map <- '--map' %in% args
    snvs <- '--snvs' %in% args
    coclus <- '--coclus' %in% args
}

snv_dir <- ''
if (snvs) {
    snv_dir <- 'snvs/'
}

require(combinat)
require(RColorBrewer)
require(gplots)
require(ggplot2)
library(grid)
require(gridExtra)
require(reshape)
require(gtools)
library(plyr)
library(gtools)

############################################################################################################################
# Helper functions
############################################################################################################################

adjust_vafs <- function(dat, pur) {
    mut_prop <- dat$prop_chrs_bearing_mutation
    vaf <- dat$adjusted_vaf
    frac <- dat$cn_frac

    cn_v <- dat$most_likely_variant_copynumber * frac
    cn_r <- dat$most_likely_ref_copynumber * (1 - frac)
    cn_r <- sapply(cn_r,function(x){if(x==1){0}else{x}})

    pn <- (1 - pur) * 2
    pr <- pur * cn_r
    pv <- pur * cn_v

    norm_const <- pn + pr + pv
    pv <- pv / norm_const

    prob <- pv * mut_prop
    adj_vaf <- (1 / prob) * vaf
    adj_vaf[adj_vaf > 2] <- 2

    return(adj_vaf)
}

get_frac <- function(x, snvs) {
    variant_cn <- x['most_likely_variant_copynumber']
    side <- x['preferred_side']
    gtype <- x['gtype']
    if (!snvs) {
        if(side==0) {
            gtype <- x['gtype1']
        } else {
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

get_runs <- function(wd) {
    runs <- c()
    for (dir in list.dirs(wd)) {
        cur_dir <- strsplit(dir,'/')[[1]]
        if (length(cur_dir) <= 1)
            next
        cur_dir <- cur_dir[length(cur_dir)]
        if (substring(cur_dir,1,3) == 'run') {
            runs <- c(runs,cur_dir)
        }
    }
    return(unique(runs))
}

get_ic_table <- function(wd, base_name, runs, allowed_ics = c('BIC', 'AIC', 'AICc', 'svc_IC'), clus_penalty = 4, snvs = FALSE) {
    snv_pref <- ''
    samp_dir <- paste(wd, base_name, sep='/')

    if (snvs) {snv_pref <- 'snvs/'}
    ic_table <- NULL
    for (run in runs) {
        ic <- read.table(paste(samp_dir, '/', run, '/', snv_pref, base_name, '_fit.txt', sep=''), sep='\t', header=F, stringsAsFactors = F)
        sc <- read.table(paste(samp_dir, '/', run, '/', base_name, '_subclonal_structure.txt', sep=''), sep='\t', header=T, stringsAsFactors = F)

        accept_solution <- !(nrow(sc) == 1 & sc[1,'CCF'] < 0.9)
        ic$V3 <- accept_solution

        ic <- ic[ic$V1%in%allowed_ics,]
        ic <- cbind(run=run, ic)
        ic_table <- rbind(ic_table, ic)
    }
    ic_table$run <- factor(ic_table$run, levels=mixedsort(unique(as.character(ic_table$run))))
    return(ic_table)
}

get_best_run <- function(wd, base_name, ic_metric, clus_penalty = 4) {
    runs <- get_runs(wd)
    ic <- get_ic_table(wd, base_name, runs, clus_penalty = clus_penalty)

    min_ic <- min(ic[ic$V1==ic_metric,'V2'])
    best_run <- as.character(ic[ic$V2==min_ic & ic$V1==ic_metric,'run'])[1]
    return(best_run)
}
get_purity <- function(wd, base_name) {
    pur <- read.table(paste(wd, base_name, '/purity_ploidy.txt', sep=''),
                      header=T, sep='\t', stringsAsFactors=F)$purity
    return(pur)
}

get_gtype <- function(x) {
    if(as.numeric(x['frac1_A'])==1) {
        return(paste(x['nMaj1_A'], x['nMin1_A'], x['frac1_A'],sep=','))
    } else {
        return(paste(paste(x['nMaj1_A'], x['nMin1_A'], x['frac1_A'],sep=','),
                     paste(as.character(x['nMaj2_A']), as.character(x['nMin2_A']), as.character(x['frac2_A']),sep=','),sep='|'))
    }
}

############################################################################################################################
# Plotting functions
############################################################################################################################

get_run_info <- function(wd, base_name, run, snvs = FALSE) {
    snv_pref <- ''
    if (snvs) {snv_pref <- 'snvs/'}
    scs_file <-  paste(wd, base_name, '/', run, '/', snv_pref, base_name, '_subclonal_structure.txt', sep = '')
    scs <- read.table(scs_file, sep = '\t', header = T)
    #     scs <- scs[scs$n_ssms>3,]
    scs <- scs[order(scs$CCF, decreasing=T), ]
    scs$new_cluster <- 1:nrow(scs)

    sv_df <- ''
    if (snvs) {
        sv_df <- read.table(paste(wd, base_name, '/', base_name, '_filtered_snvs.tsv', sep=''), header=T, sep='\t', stringsAsFactors=F)
    } else {
        sv_df <- read.table(paste(wd, base_name, '/', base_name, '_filtered_svs.tsv', sep=''), header=T, sep='\t', stringsAsFactors=F)
    }
    pur <- read.table(paste(wd, base_name, '/purity_ploidy.txt', sep=''),
                      header=T, sep='\t', stringsAsFactors=F)$purity

    cc_file <-  paste(wd, base_name, '/', run, '/', snv_pref, base_name, '_cluster_certainty.txt', sep = '')
    cc <- read.table(cc_file, sep = '\t', stringsAsFactors = F, header = T)

    mlcn_file <- paste(wd, base_name, '/', run, '/', snv_pref, base_name, '_most_likely_copynumbers.txt', sep='')
    mlcn <- read.table(mlcn_file, header = T, sep = '\t', stringsAsFactors = F)

    merge_cols <- c('chr1', 'pos1', 'dir1', 'chr2', 'pos2', 'dir2')
    if (snvs) {
        colnames(sv_df)[1] <- 'chr'
        colnames(cc)[1] <- 'chr'
        merge_cols <- c('chr', 'pos')
    }
    dat <- merge(sv_df, mlcn, by.x=merge_cols, by.y=merge_cols)
    dat <- merge(dat, cc, by.x=merge_cols, by.y=merge_cols)
    dat$cn_frac <- apply(dat, 1, function(x){get_frac(x, snvs)})
    if(snvs) {dat$adjusted_vaf <- dat$var / (dat$ref + dat$var)}
    dat <- cbind(dat, CCF=adjust_vafs(dat, pur))

    dat$cluster <- NA
    for (j in 1:nrow(scs)) {
        dat[dat$most_likely_assignment==scs$cluster[j], 'cluster'] <- scs$new_cluster[j]
    }
    if (snvs) {
        dat$sv <- paste(dat$chr, dat$pos, sep='_')
    } else {
        dat$sv <- paste(dat$chr1, dat$pos1, dat$dir1, dat$chr2, dat$pos2, dat$dir2, sep=':')
    }
    dat <- dat[!duplicated(dat),]
    scs$cluster <- scs$new_cluster
    scs <- scs[,1:4]
    scs$variant_proportion <- scs$n_ssms/sum(scs$n_ssms)
    return(list(dat,scs))
}

plot_ccf_hist <- function(wd, base_name, snvs, pick_run='best') {
    if (pick_run=='best') {
        pick_run <- get_best_run(wd, base_name, 'svc_IC')
    }

    x <- get_run_info(wd, base_name, pick_run, snvs)
    dat <- x[[1]]
    sc <- x[[2]]
    sc$cluster <- sc$cluster-1
    sc$clus_renum <- rank(sc$CCF)

    dat$clus_renum <- NULL
    for (i in 1:nrow(sc)) {
        dat[dat$most_likely_assignment==sc$cluster[i], 'clus_renum'] <- sc[i,'clus_renum']
    }

    pp <- read.delim(paste(wd, base_name, '/purity_ploidy.txt', sep=''))
    pur <- pp$purity

    above_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) >= 0.04
    below_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) < 0.04
    clus_intercepts <- 1 / pur * as.numeric(sc$proportion[above_ssm_th & sc$n_ssms > 2])
    clus_intercepts_minor <- 1 / pur * as.numeric(sc$proportion[below_ssm_th | sc$n_ssms<=2])

    clus_freqs <- table(dat$most_likely_assignment)
    minor_clusts <- names(clus_freqs[clus_freqs/nrow(dat) <= 0.04])

    if (length(minor_clusts) != 0) {
        dat <- dat[!dat$most_likely_assignment%in%minor_clusts,]
    }

    ccf_hist <- ggplot(dat, aes(x=as.numeric(dat$CCF),
                                fill=factor(most_likely_assignment),color=factor(most_likely_assignment))) +
        theme_minimal() +
        xlim(0,2) + geom_histogram(alpha=0.3,position='identity',binwidth=0.05)+xlab('CCF') +
        geom_vline(xintercept=clus_intercepts, colour='blue', size=1) + ylab('') +
        scale_fill_brewer(palette = 'Set1', name = "Cluster") +
        scale_color_brewer(palette = 'Set1', name = "Cluster") +
        theme(plot.title = element_text(size = 16, face = "bold"),
              axis.title = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14)) +
        ggtitle(pick_run)

    if (length(clus_intercepts_minor) > 0) {
        ccf_hist <- ccf_hist + geom_vline(xintercept=clus_intercepts_minor,colour='red',lty=2)
    }
    return(ccf_hist)
}




# Get number of runs
runs <- get_runs('.')
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

# get_adjust_factor <- function(dat, pur) {
#     cn_v <- dat$most_likely_variant_copynumber
#     cn_r <- dat$most_likely_ref_copynumber
#     mut_prop <- dat$prop_chrs_bearing_mutation
#     vaf <- dat$adjusted_vaf
#     frac <- dat$cn_frac
#     # prob <- (cn_v * mut_prop * pur) / (2 * (1 - pur) + cn_v * pur)
#     prob <- (cn_v * mut_prop * pur) / (2 * (1 - pur) + (cn_v * pur * frac) + (cn_r * pur * (1. - frac)))
#     return((1 / prob))
# }

# gg_color_hue <- function(n) {
#     hues = seq(15, 375, length=n+1)
#     hcl(h=hues, l=65, c=100)[1:n]
# }

# ggQQ <- function(dat) {
#     p <- ggplot(dat) +
#         stat_qq(aes(sample=CCF, colour = factor(most_likely_assignment)), alpha = 0.5)
#
#     dat_tmp <- dat[!is.na(dat$CCF),]
#
#     clusts <- as.numeric(names(table(dat_tmp$most_likely_assignment)))
#     cols <- gg_color_hue(length(clusts))
#
#     for (i in 1:length(clusts)) {
#         clus <- clusts[i]
#         tmp <- dat_tmp[dat_tmp$most_likely_assignment == clus, 'CCF']
#         y <- quantile(tmp, c(0.25, 0.75))
#         x <- qnorm(c(0.25, 0.75))
#         slope <- diff(y)/diff(x)
#         intercept <- y[1L] - slope * x[1L]
#
#         p <- p + geom_abline(slope = slope, intercept = intercept, color=cols[i], alpha=0.5)
#     }
#
#     return(p)
# }

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
    plot1 <- plot_ccf_hist('../', id, snvs, run)

    x <- get_run_info('../', id, run, snvs = snvs)
    dat <- x[[1]]
    plot2 <- ggplot(dat, aes(x=CCF, y=adjusted_vaf,
                             fill=factor(most_likely_assignment),
                             colour=factor(most_likely_assignment))) +
                    scale_fill_brewer(palette = 'Set1', name = "Cluster") +
                    scale_color_brewer(palette = 'Set1', name = "Cluster") +
                    theme(plot.title = element_text(size = 16, face = "bold"),
                          axis.title = element_text(size = 16),
                          axis.text.x = element_text(size = 14),
                          axis.text.y = element_text(size = 14)) +
                    geom_point(size=1) + xlim(0,2) + ylim(0,1) + ylab('VAF')

    # attach table for convenience, also add BIC/AIC
    sv_clust <- read.table(paste(run, '/', snv_dir, id, '_subclonal_structure.txt', sep=''),
                           header=T, sep='\t', stringsAsFactors=F)
    tabout <- sv_clust[order(as.numeric(sv_clust[,3]),decreasing=TRUE),]
    sc_tab <- tableGrob(tabout, rows=c())

    ic_tab <- NULL
    if (map) {
        ic <- ic_table[ic_table$run==run,]
        rownames(ic) <- ic$V1
        ic <- ic[,-c(1:2)]
        colnames(ic)[1] <- 'value'
        ic$value <- round(ic$value, 4)
        ic_tab <- tableGrob(ic)
    }

    if (coclus) {
        plot3 <- plot_ccf_hist('../', id, snvs = TRUE, pick_run = run)
        dat <- get_run_info('../', id, run, snvs = TRUE)[[1]]
        plot4 <- ggplot(dat, aes(x=CCF, y=adjusted_vaf,
                                 fill=factor(most_likely_assignment),
                                 colour=factor(most_likely_assignment))) +
            scale_fill_brewer(palette = 'Set1', name = "Cluster") +
            scale_color_brewer(palette = 'Set1', name = "Cluster") +
            theme(plot.title = element_text(size = 16, face = "bold"),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14)) +
            geom_point(size=1) + xlim(0,2) + ylim(0,1) + ylab('VAF')

        if (map) {
            height <- 8+round(nrow(tabout)*0.2)
            pdf(paste(id, run, 'fit.pdf',sep='_'), height=height, width=9)
            grid.arrange(sc_tab, ic_tab,
                         plot1 + ggtitle(paste(run, 'SVs')), plot2 + ggtitle('SVs'),
                         plot3 + ggtitle(paste(run, 'SNVs')), plot4 + ggtitle('SNVs'), nrow=3)
            dev.off()
        } else {
            blank <- grid.rect(gp=gpar(col="white"))
            height <- 8+round(nrow(tabout)*0.2)
            pdf(paste(id, run, 'fit.pdf',sep='_'), height=height, width=9)
            grid.arrange(sc_tab, blank, plot1 + ggtitle(paste(run, 'SVs')), plot2 + ggtitle('SVs'),
                         plot3 + ggtitle(paste(run, 'SNVs')), plot4 + ggtitle('SNVs'), nrow=3)
            dev.off()
        }
    } else {
        if (map) {
            height <- 5+round(nrow(tabout)*0.2)
            pdf(paste(id, run, 'fit.pdf',sep='_'), height=height, width=9)
            grid.arrange(sc_tab, plot1, ic_tab, plot2, ncol=2)
            dev.off()
        } else {
            height <- 7+round(nrow(tabout)*0.2)
            pdf(paste(id, run, 'fit.pdf',sep='_'), height=height, width=5)
            grid.arrange(sc_tab, plot1, plot2, ncol=1)
            dev.off()
        }
    }
}

var_id <- c('chr1','pos1','chr2','pos2')
if (snvs) {var_id <- c('chr', 'pos')}

all_runs_ccfs <- NULL
all_runs_scs <- NULL
for(run in runs) {
    runinf <- get_run_info('../', id, run, snvs)
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
    pdf(pdf_name, height = 6, width = 18 + (length(runs)/2))
    grid.arrange(rp, ic_plot, svic_plot, nrow = 1)
    dev.off()
} else {
    pdf_name <- paste(id, '_run_summary.pdf', sep='')
    pdf(pdf_name, height = 6, width = 8 + (length(runs)/6))
    print(rp)
    dev.off()
}
