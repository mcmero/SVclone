############################################################################################################################
# SVclone output post processing
#
# Draw some helpful plots to interpret SVclone output
# plots run summary, IC and SNV/SV histograms
############################################################################################################################

args <- commandArgs(trailingOnly = TRUE)

if(length(args)==0 | args[1]=='-h') {
    print('Usage: Rscript <workingdir> <sample id> <battenberg CNs> (opt) --map --snvs --coclus')
}

rundir <- getwd()
wd <- args[1]
setwd(wd)

id <- args[2]
bbf <- args[3]
map <- FALSE
clus_stab <- FALSE
snvs <- FALSE
coclus <- FALSE
allowed_ics <- c('svc_IC')
snv_dir <- ''

map <- '--map' %in% args
snvs <- '--snvs' %in% args
coclus <- '--coclus' %in% args

if (snvs) {snv_dir <- 'snvs/'}
if (coclus) {snvs<-FALSE}

suppressMessages(library(combinat))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(reshape))
suppressMessages(library(gtools))
suppressMessages(library(plyr))
suppressMessages(library(gtools))

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
    if (gtype==''){return(1)}
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

get_runs <- function(wd, filter_out_post=TRUE) {
    runs <- c()
    dirlist <- list.dirs(wd)
    pa <- dirlist[grep('post_assign', dirlist)]
    if (length(pa)>0) {dirlist <- dirlist[!dirlist%in%pa]}

    for (dir in dirlist) {
        cur_dir <- strsplit(dir,'/')[[1]]
        if (length(cur_dir) <= 1)
            next
        cur_dir <- cur_dir[length(cur_dir)]
        if (substring(cur_dir,1,3) == 'run') {
            runs <- c(runs,cur_dir)
        }
    }
    if (!filter_out_post & length(pa) > 0) {
        pa <- pa[sapply(pa, function(x){length(strsplit(x, './')[[1]])<3})]
        pa <- as.character(sapply(pa, function(x){x<-strsplit(x,'/')[[1]];x<-x[length(x)]}))
        runs <- c(runs, pa)
    }
    return(unique(runs))
}

get_ic_table <- function(wd, base_name, runs, allowed_ics = c('BIC', 'AIC', 'AICc', 'svc_IC'), clus_penalty = 4, snvs = FALSE) {
    snv_pref <- ''
    samp_dir <- wd

    if (snvs) {snv_pref <- 'snvs/'}
    ic_table <- NULL
    for (run in runs) {
        ic <- read.table(paste(samp_dir, '/', run, '/', snv_pref, base_name, '_fit.txt', sep=''), sep='\t', header=F, stringsAsFactors = F)
        sc <- read.table(paste(samp_dir, '/', run, '/', snv_pref, base_name, '_subclonal_structure.txt', sep=''), sep='\t', header=T, stringsAsFactors = F)

        accept_solution <- !(nrow(sc) == 1 & sc[1,'CCF'] < 0.9)
        ic$V3 <- accept_solution

        ic <- ic[ic$V1%in%allowed_ics,]
        ic <- cbind(run=run, ic)
        ic_table <- rbind(ic_table, ic)
    }
    ic_table$run <- factor(ic_table$run, levels=mixedsort(unique(as.character(ic_table$run))))
    return(ic_table)
}

get_best_run <- function(wd, base_name, ic_metric, clus_penalty = 4, snvs = FALSE) {
    runs <- get_runs(wd)
    ic <- get_ic_table(wd, base_name, runs, clus_penalty = clus_penalty, snvs = snvs)

    min_ic <- min(ic[ic$V1==ic_metric,'V2'])
    best_run <- as.character(ic[ic$V2==min_ic & ic$V1==ic_metric,'run'])[1]
    return(best_run)
}
get_purity <- function(wd, base_name) {
    pur <- read.table(paste(wd, '/purity_ploidy.txt', sep=''),
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

get_run_info <- function(wd, base_name, run, snvs = FALSE, post = FALSE) {
    snv_pref <- ''
    if (snvs) {snv_pref <- 'snvs/'}
    scs_file <-  paste(wd, '/', run, '/', snv_pref, base_name, '_subclonal_structure.txt', sep = '')
    scs <- read.table(scs_file, sep = '\t', header = T)
    colnames(scs)[2] <- 'n_ssms'
    scs <- scs[,c('cluster', 'n_ssms', 'proportion', 'CCF')]
    scs <- scs[order(scs$CCF, decreasing=T), ]
    scs$new_cluster <- 1:nrow(scs)

    sv_df <- ''
    if (snvs) {
        postfix <- if(!post){c('_filtered_snvs.tsv')} else {c('_filtered_snvs_post_assign.tsv')}
        sv_df <- read.table(paste(wd, '/', base_name, postfix, sep=''), header=T, sep='\t', stringsAsFactors=F)
    } else {
        postfix <- if(!post){c('_filtered_svs.tsv')} else {c('_filtered_svs_post_assign.tsv')}
        sv_df <- read.table(paste(wd, '/', base_name, postfix, sep=''), header=T, sep='\t', stringsAsFactors=F)
    }
    pur <- get_purity(wd, base_name)

    cc_file <-  paste(wd, '/', run, '/', snv_pref, base_name, '_cluster_certainty.txt', sep = '')
    cc <- read.table(cc_file, sep = '\t', stringsAsFactors = F, header = T)

    mlcn_file <- paste(wd, '/', run, '/', snv_pref, base_name, '_most_likely_copynumbers.txt', sep='')
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

plot_hist <- function(wd, base_name, snvs, pick_run='best', varclass=FALSE, vaf=FALSE, post=FALSE) {
    if (pick_run=='best' & map) {
        pick_run <- best_run
    }

    x <- get_run_info(wd, base_name, pick_run, snvs, post)
    dat <- x[[1]]
    sc <- x[[2]]
    sc$cluster <- sc$cluster-1
    sc$clus_renum <- rank(sc$CCF)

    dat$clus_renum <- NULL
    for (i in 1:nrow(sc)) {
        dat[dat$most_likely_assignment==sc$cluster[i], 'clus_renum'] <- sc[i,'clus_renum']
    }
    pur <- get_purity(wd, base_name)

    above_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) >= 0.04
    below_ssm_th <- sc$n_ssms / (sum(sc$n_ssms)) < 0.04
    clus_intercepts <- 1 / pur * as.numeric(sc$proportion[above_ssm_th & sc$n_ssms > 2])
    clus_intercepts_minor <- 1 / pur * as.numeric(sc$proportion[below_ssm_th | sc$n_ssms<=2])

    clus_freqs <- table(dat$most_likely_assignment)
    minor_clusts <- names(clus_freqs[clus_freqs/nrow(dat) <= 0.04])

    # if (length(minor_clusts) != 0) {
    #     dat <- dat[!dat$most_likely_assignment%in%minor_clusts,]
    # }

    dat$most_likely_assignment <- factor(dat$most_likely_assignment)
    plotvar <- 'Cluster'
    if (varclass) {
        plotvar <- 'Classification'
        if (!'classification'%in%colnames(dat)) {
            dat$classification <- rep('SNV', nrow(dat))
        }
        if (vaf) {
            var_hist <- ggplot(dat, aes(x=adjusted_vaf,fill=classification,color=classification)) +
                geom_histogram(alpha=0.3,binwidth=0.05) + xlab('VAF')
        } else {
            var_hist <- ggplot(dat, aes(x=CCF,fill=classification,color=classification)) +
                geom_histogram(alpha=0.3,binwidth=0.05) + xlab('CCF')
        }
    } else {
        if (vaf) {
            var_hist <- ggplot(dat, aes(x=adjusted_vaf,fill=most_likely_assignment,color=most_likely_assignment)) +
                geom_histogram(alpha=0.3,position='identity',binwidth=0.05) + xlab('VAF')
        } else {
            var_hist <- ggplot(dat, aes(x=CCF,fill=most_likely_assignment,color=most_likely_assignment)) +
                geom_histogram(alpha=0.3,position='identity',binwidth=0.05) + xlab('CCF')
        }
    }

    var_hist <- var_hist + theme_minimal() + xlim(0,2) + ggtitle(pick_run) + ylab('') +
        scale_fill_brewer(palette = 'Set1', name = plotvar) +
        scale_color_brewer(palette = 'Set1', name = plotvar) +
        theme(plot.title = element_text(size = 16, face = 'bold'),
              axis.title = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14))

    if (!vaf) {
        var_hist <- var_hist + geom_vline(xintercept=clus_intercepts, colour='blue', size=1)
        if (length(clus_intercepts_minor) > 0) {
            var_hist <- var_hist + geom_vline(xintercept=clus_intercepts_minor,colour='red',lty=2)
        }
    }
    return(var_hist)
}

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

############################################################################################
# Cluster proportions
############################################################################################

runs <- get_runs('./')
clusts <- NULL
for (run in runs) {
    sv_clust <- read.table(paste(run, '/', snv_dir, id, '_subclonal_structure.txt', sep=''),
                           header=T, sep='\t', stringsAsFactors=F)
    clusts <- rbind(clusts, data.frame(run=run, n_ssms=sv_clust[,2], cluster=sv_clust$cluster))
}

pdf(paste(id, '_cluster_hist.pdf',sep=''),height=4, width=max(3,length(runs)*0.6))
ggplot(clusts, aes(x=factor(run), y=n_ssms, fill=factor(cluster))) + geom_bar(stat='identity') + theme_minimal()
dev.off()

############################################################################################
# Plot AIC and BIC for runs
############################################################################################

if (length(args)>2 & map) {
    #print('Plotting AIC & BIC metrics...')
    ic_table <- NULL
    for (run in runs) {
        ic <- read.table(paste(run, '/', snv_dir, id, '_fit.txt', sep=''), sep='\t', header=F, stringsAsFactors = F)
        ic <- ic[ic$V1%in%allowed_ics,]
        ic <- cbind(run=run, ic)
        ic_table <- rbind(ic_table, ic)
    }

    colnames(ic_table)[2:3] <- c('IC', 'value')
    ic_table$run <- factor(ic_table$run, levels=mixedsort(unique(as.character(ic_table$run))))
    ic_table$IC <- factor(ic_table$IC, levels=allowed_ics)

#    g1 <- ggplot(ic_table[ic_table$IC=='BIC',],
#                 aes(y=value, x=run, group=IC)) + ylab('value') + geom_line() + ggtitle('BIC') +
#                 theme_minimal()
#    g2 <- ggplot(ic_table[ic_table$IC=='AIC',],
#                 aes(y=value, x=run, group=IC)) + ylab('value') + geom_line() + ggtitle('AIC') +
#                 theme_minimal()
#    g3 <- ggplot(ic_table[ic_table$IC=='AICc',],
#                 aes(y=value, x=run, group=IC)) + ylab('value') + geom_line() + ggtitle('AICc') +
#                 theme_minimal()
#    g4 <- ggplot(ic_table[ic_table$IC=='svc_IC',],
#                 aes(y=value, x=run, group=IC)) + ylab('value') + geom_line() + ggtitle('svc_IC') +
#                 theme_minimal()
#    pdf(paste(id, 'ic_plots.pdf',sep='_'), height=6)
#    grid.arrange(g1,g2,g3,g4)
#    dev.off()

    ic_table <- cast(ic_table, run~IC)

#    min_bic <- ic_table[min(ic_table$BIC)==ic_table$BIC,]
#    min_bic$AIC <- min_bic$run
#    min_bic$run <- 'min_BIC'
#    min_bic$AICc <- min_bic$BIC
#    min_bic$BIC <- NA; min_bic$svc_IC <- NA
#
#    min_aic <- ic_table[min(ic_table$AIC)==ic_table$AIC,]
#    min_aic$AICc <- min_aic$AIC
#    min_aic$AIC <- min_aic$run
#    min_aic$run <- 'min_AIC'
#    min_aic$BIC <- NA; min_aic$svc_IC <- NA
#
#    min_aicc <- ic_table[min(ic_table$AICc)==ic_table$AICc,]
#    min_aicc$AIC <- min_aicc$run
#    min_aicc$run <- 'min_AICc'
#    min_aicc$BIC <- NA; min_aicc$svc_IC <- NA

    min_svcic <- ic_table[min(ic_table$svc_IC,na.rm=T)==ic_table$svc_IC,]
#    min_svcic$AIC <- min_svcic$run
    min_svcic$run <- 'min_svc_IC'
#    min_svcic$AICc <- min_svcic$svc_IC
#    min_svcic$BIC <- NA; min_svcic$svc_IC <- NA

#    ic_table <- rbind(ic_table, min_bic)
#    ic_table <- rbind(ic_table, min_aic)
#    ic_table <- rbind(ic_table, min_aicc)
    ic_table <- rbind(ic_table, min_svcic)

    write.table(ic_table, paste(id,'_ic_metrics.csv', sep=''), sep=',', quote=F, row.names=F, na = '')
}

############################################################################################
# Plot histogram of clusters
############################################################################################

if (map) {
    ic_table <- NULL
    for (run in runs) {
        ic <- read.table(paste(run, '/', snv_dir, id, '_fit.txt', sep=''),
                         sep='\t', header=F, stringsAsFactors = F)
        ic <- ic[ic$V1%in%allowed_ics,]
        ic <- cbind(run=run, ic)
        ic_table <- rbind(ic_table, ic)
    }
    ic_table$rank <- NULL
    for (aic in allowed_ics) {
        ic_table[ic_table$V1 == aic, 'rank'] <-
            rank(ic_table[ic_table$V1 == aic, 'V2'])
    }
    colnames(ic_table)[2:3] <- c('IC', 'value')
    ic_table$run <- factor(ic_table$run, levels=mixedsort(unique(as.character(ic_table$run))))
    ic_table$IC <- factor(ic_table$IC, levels=allowed_ics)
}

best_run <- 'none'
if(map) {
    best_run <- get_best_run('./', id, 'svc_IC', snvs = snvs)
    post_assign1 <- 'best_run_snvs_post_assign'
    post_assign2 <- 'best_run_svs_post_assign'
    post_assign3 <- 'best_run_coclus_post_assign'
    if (file.exists(post_assign1)){best_run <- post_assign1}
    if (file.exists(post_assign2)){best_run <- post_assign2}
    if (file.exists(post_assign3)){best_run <- post_assign3}
}

runs <- get_runs('./', filter_out_post=F)
for (run in runs) {
    outname <- paste(id, run, 'fit.pdf',sep='_')
    if (map & run==best_run) {
        outname <- paste(id, run, 'best_fit.pdf',sep='_')
    }

    post <- grepl('post_assign', run)
    plot1 <- plot_hist('./', id, snvs, run, post = post)
    plot2 <- plot_hist('./', id, snvs = snvs, pick_run = run, vaf = TRUE, post = post)
    x <- get_run_info('./', id, run, snvs = snvs, post = post)
    dat <- x[[1]]

    plot3 <- ggplot(dat, aes(x=CCF, y=adjusted_vaf,
                             fill=factor(most_likely_assignment),
                             colour=factor(most_likely_assignment))) +
                    scale_fill_brewer(palette = 'Set1', name = 'Cluster') +
                    scale_color_brewer(palette = 'Set1', name = 'Cluster') +
                    theme(plot.title = element_text(size = 16, face = 'bold'),
                          axis.title = element_text(size = 16),
                          axis.text.x = element_text(size = 14),
                          axis.text.y = element_text(size = 14)) +
                    geom_point(size=1) + xlim(0,2) + ylim(0,1) + ylab('VAF') + theme_minimal()

    # attach table for convenience, also add BIC/AIC
    sv_clust <- read.table(paste(run, '/', snv_dir, id, '_subclonal_structure.txt', sep=''),
                           header=T, sep='\t', stringsAsFactors=F)
    colnames(sv_clust)[2] <- 'n_ssms'

    ic_tab <- NULL
    if (map & length(grep('post', run))==0) {
        ic <- ic_table[ic_table$run==run,]
        rownames(ic) <- ic$IC
        ic <- ic[,-c(1:2)]
        colnames(ic)[1] <- 'value'
        ic$value <- round(ic$value, 4)
        ic_tab <- tableGrob(ic)
    }

    if (coclus) {
        snv <- get_run_info('./', id, run, snvs = TRUE, post = post)[[1]]
        svs <- get_run_info('./', id, run, snvs = FALSE, post = post)[[1]]
        svs <- data.frame(table(svs$most_likely_assignment))
        colnames(svs) <- c('cluster', 'SVs')
        sv_clust <- merge(sv_clust, svs, by='cluster', all.x=T, all.y=T)
        sv_clust[is.na(sv_clust)] <- 0
        sv_clust$SNVs <- sv_clust$n_ssms-sv_clust$SVs
        colnames(sv_clust)[2] <- 'variants'
        tabout <- sv_clust[order(as.numeric(sv_clust[,3]),decreasing=TRUE),]
        tabout <- tabout[,c('cluster', 'SNVs', 'SVs', 'proportion', 'CCF')]
        sc_tab <- tableGrob(tabout, rows=c())

        plot4 <- ggplot(snv, aes(x=CCF, y=adjusted_vaf,
                                 fill=factor(most_likely_assignment),
                                 colour=factor(most_likely_assignment))) +
            scale_fill_brewer(palette = 'Set1', name = 'Cluster') +
            scale_color_brewer(palette = 'Set1', name = 'Cluster') +
            theme(plot.title = element_text(size = 16, face = 'bold'),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14)) +
            geom_point(size=1) + xlim(0,2) + ylim(0,1) + ylab('VAF') + theme_minimal()
        plot5 <- plot_hist('./', id, snvs = TRUE, pick_run = run, post = post)
        plot6 <- plot_hist('./', id, snvs = FALSE, pick_run = run, vaf = TRUE, post = post)
        plot7 <- plot_hist('./', id, snvs = TRUE, pick_run = run, varclass = TRUE, post = post)
        plot8 <- plot_hist('./', id, snvs = FALSE, pick_run = run, varclass = TRUE, post = post)

        if (map & length(grep('post', run))==0) {
            height <- 12+round(nrow(tabout)*0.2)
            pdf(outname, height=height, width=9)
            grid.arrange(sc_tab, ic_tab,
                         plot1 + ggtitle(paste(run, 'SVs')), plot5 + ggtitle('SNVs'),
                         plot2 + ggtitle(paste(run, 'SVs')), plot6 + ggtitle('SNVs'),
                         plot8 + ggtitle(paste(run, 'SVs')), plot7 + ggtitle(paste(run, 'SNVs')),
                         plot3 + ggtitle(paste(run, 'SVs')), plot4 + ggtitle(paste(run, 'SNVs')), ncol=2)
            dev.off()
        } else {
            height <- 12+round(nrow(tabout)*0.2)
            pdf(outname, height=height, width=9)
            grid.arrange(sc_tab, grid.rect(gp=gpar(col='white')),
                         plot1 + ggtitle(paste(run, 'SVs')), plot5 + ggtitle('SNVs'),
                         plot2 + ggtitle(paste(run, 'SVs')), plot6 + ggtitle('SNVs'),
                         plot8 + ggtitle(paste(run, 'SVs')), plot7 + ggtitle(paste(run, 'SNVs')),
                         plot3 + ggtitle(paste(run, 'SVs')), plot4 + ggtitle(paste(run, 'SNVs')), ncol=2)
            dev.off()
        }
    } else {
        varname <- 'SVs'
        if(snvs){varname<-'SNVs'}
        colnames(sv_clust)[2] <- varname
        tabout <- sv_clust[order(as.numeric(sv_clust[,3]),decreasing=TRUE),]
        tabout <- tabout[,c('cluster', varname, 'proportion', 'CCF')]
        sc_tab <- tableGrob(tabout, rows=c())
        plot4 <- plot_hist('./', id, snvs, run, varclass = TRUE, post = post)
        if (map & !post) {
            height <- 7+round(nrow(tabout)*0.2)
            pdf(outname, height=height, width=9)
            grid.arrange(sc_tab, ic_tab,
                         plot1, plot2,
                         plot4, plot3, ncol=2)
            dev.off()
        } else {
            height <- 10+round(nrow(tabout)*0.2)
            pdf(outname, height=height, width=4)
            grid.arrange(sc_tab, plot1, plot2, plot4, plot3, ncol=1)
            dev.off()
        }
    }
}

############################################################################################
# Summary plots
############################################################################################

var_id <- c('chr1','pos1','chr2','pos2')
if (snvs) {var_id <- c('chr', 'pos')}

all_runs_ccfs <- NULL
all_runs_scs <- NULL
runs <- get_runs('./')
for(run in runs) {
    runinf <- get_run_info('./', id, run, snvs)
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
    scale_size(range = c(4,25), guide = FALSE) + ggtitle('Clustering summary') +
    geom_point(aes(x = run, y = CCF, size=variant_proportion), colour = '#4d4d4d', alpha=0.8) +
    theme_minimal()

if (map) {
    svic_plot <- ggplot(ic_table[ic_table$IC=='svc_IC',], aes(y=value, x=run, group=1)) +
        ylab('value') + geom_line() + ggtitle('SVclone IC') + theme_minimal()
    pdf_name <- paste(id, '_run_summary.pdf', sep='')
    pdf(pdf_name, height = 6, width = 14 + (length(runs)/2))
    suppressWarnings(grid.arrange(rp, svic_plot, nrow = 1))
    dev.off()
} else {
    pdf_name <- paste(id, '_run_summary.pdf', sep='')
    pdf(pdf_name, height = 6, width = 7 + (length(runs)/6))
    suppressWarnings(print(rp))
    dev.off()
}

############################################################################################
# Circos plots
############################################################################################

setwd(rundir)
if (!grepl('^--', bbf) & file.exists(bbf)) {
    print('Plotting circos...')
    bb <- read.delim(bbf, sep='\t', stringsAsFactors = F)

    if (length(grep('battenberg', colnames(bb))) > 0) {
        colnames(bb)[2:3] <- c('startpos', 'endpos')
        colnames(bb)[grep('nMaj1_A', colnames(bb))] <- 'nMaj1_A'
        colnames(bb)[grep('nMin1_A', colnames(bb))] <- 'nMin1_A'
        colnames(bb)[grep('nMaj2_A', colnames(bb))] <- 'nMaj2_A'
        colnames(bb)[grep('nMin2_A', colnames(bb))] <- 'nMin2_A'
    } else if (length(grep('chromosome', colnames(bb))) > 0) {
        # consensus format
        bb <- bb[,1:8]
        colnames(bb) <- c('chr', 'startpos', 'endpos', 'total_cn', 'nMaj1_A', 'nMin1_A', 'star', 'level')
        bb$nMaj2_A <- NA; bb$nMin2_A <- NA
    } else if (ncol(bb)==1) {
        # ASCAT caveman format
        bb <- read.delim(bbf, sep=',', stringsAsFactors = F, header = F)
        colnames(bb) <- c('ID','chr','startpos','endpos','norm_total','norm_minor','total_cn','nMin1_A')
        bb$nMaj1_A <- bb$total_cn - bb$nMin1_A
        bb$nMaj2_A <- NA; bb$nMin2_A <- NA
    }

    clon <- bb[is.na(bb$nMaj2_A),]
    subclon <- bb[!is.na(bb$nMaj2_A),]

    pdat <- clon[,c('chr','startpos','endpos')]
    pdat <- cbind(pdat, value=clon$nMaj1_A+clon$nMin1_A)
    pdat$chr <- paste('chr', pdat$chr,sep='')
    colnames(pdat) <- c('chr','start','end','value')
    pdat$value[pdat$value > 6] <- 6

    if (nrow(subclon) > 0) {
        pdat2 <- subclon[,c('chr','startpos','endpos')]
        pdat2 <- cbind(pdat2, value=apply(cbind(subclon$nMaj2_A+subclon$nMin2_A,subclon$nMaj1_A+subclon$nMin1_A),1,mean))
        pdat2$chr <- paste('chr',pdat2$chr,sep='')

        pdat2$value[pdat2$value > 6] <- 6
        colnames(pdat2)<-c('chr','start','end','value')

        pdat <- list(pdat, pdat2)
    }
    colours <- c('#0000FF80','#FF000080','darkgreen','#0000FF40','#FF000040','#00FF0040')

    run <- 'run0'; if (map) {run <- best_run}
    tmp <- get_run_info(wd, id, run, snvs)
    dat <- tmp[[1]]
    sv_clust <- tmp[[2]]
    sv_clust$cluster <- sv_clust$cluster-1

    suppressMessages(require(circlize))
    pdf(paste(wd, '/', id, '_', run, '_circos.pdf', sep=''), height=12, width=12)
    par(mar = c(1, 1, 1, 1))
    circos.initializeWithIdeogram(plotType = c('axis','labels'))

    if (coclus | snvs) {
        res_snvs <- read.table(paste(wd, '/', run, '/snvs/', id, '_cluster_certainty.txt', sep=''),
                               sep='\t', header=T, stringsAsFactors=F)
        psnvs <- list()
        count <- 1
        for(i in 1:nrow(sv_clust)) {
            curr <- res_snvs[res_snvs$most_likely_assignment==sv_clust[i,1],]
            if(nrow(curr)>1)
            {
                bed <- curr[,c(1,2)]
                bed[,1] <- paste('chr',bed[,1],sep='')
                bed <- cbind(bed,bed[,2]+1,as.numeric(sv_clust[i,3]))
                colnames(bed) <- c('chr','start','end','value')
                psnvs[[count]] <- bed
                count <- count + 1
            }
        }
        circos.genomicDensity(psnvs, col=colours, overlap=FALSE)
    }
    circos.genomicTrackPlotRegion(pdat,ylim=c(0,6),
                                  panel.fun=function(region,value,...){
                                      i=getI(...)
                                      circos.genomicLines(region,value,type='segment',lwd=3,col=colours[i],...)
                                  })
    if (coclus | !snvs) {
        for(j in 1:nrow(dat))
        {
            x <- dat[j,]
            ccf <- x$CCF
            if(ccf<0.9) {
                lcol=colours[2]
            } else {
                lcol=colours[1]
            }
            circos.link(paste('chr',as.character(x[1]),sep=''),
                        as.numeric(x[2]),
                        paste('chr',as.character(x[4]),sep=''),
                        as.numeric(x[5]),col=lcol,lwd=2)
        }
    }
    dev.off()
}
print('Finished plotting diagnostic plots!')
