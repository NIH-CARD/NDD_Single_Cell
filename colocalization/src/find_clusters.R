#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(ggplot2)
library(ggthemes)

#### ARGS ##########################################################################################
args <- commandArgs(trailingOnly=TRUE)
GWAS_FILENAME <- args[1]
RECOMBINATION_BED_FILE <- args[2]
OUTDIR=paste0(dirname(GWAS_FILENAME), '/clusters')

#### FUNCTIONS #####################################################################################

import_recombination_bed <- function(filename) {
    # Read in .bed file containing recombination rates
    # Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
    bed <- fread(file = RECOMBINATION_BED_FILE,
                header = TRUE,
                showProgress = FALSE,
                col.names = c("CHR","BP","c","cM")
    )
    bed[, c := NULL]
    setkey(bed, CHR, BP)
    bed_extension <- bed[, list('BP'=max(BP)*1.1, 'cM'=max(cM)), by=CHR]    # extend ends of chr by 10%, with same cM
    bed <- unique(rbindlist(list(bed, bed_extension)))
    setkey(bed, CHR, BP)
    return(bed)
}


build_pos_to_cM <- function(recombination_bed_table) {
# Generate functions to translate base pair position to units of cM
    #cM_to_POS <- new.env()
    env1 <- new.env()
    for(chr in unique(bed[,CHR])) {
        #cM_to_POS[[as.character(chr)]] <- approxfun(y=c(0, bed[CHR==chr][,BP]), x=c(0,bed[CHR==chr][,cM]), ties='ordered')
        env1[[as.character(chr)]] <- approxfun(x=c(0, bed[CHR==chr][,BP]), y=c(0,bed[CHR==chr][,cM]), ties='ordered')
    }
    return(env1)
}


build_cM_to_pos <- function(recombination_bed_table) {
# Generate function to translate cM to base pair position
    env1 <- new.env()
    for(chr in unique(bed[,CHR])) {
        env1[[as.character(chr)]] <- approxfun(y=c(0, bed[CHR==chr][,BP]), x=c(0,bed[CHR==chr][,cM]), ties='ordered')
    }
    return(env1)
}



get_gwas_filestem <- function(filename) {
    return(strsplit(basename(filename), split='\\.')[[1]][1])
}


get_gwas_hits <- function(gwas_filename) {
    dat <- fread(gwas_filename)
    setnames(dat, toupper(colnames(dat)))
    setkey(dat, CHR, BP)


    desired_cols <- c('SNP','CHR','BP','A1','A2','BETA','FRQ','SE','P')
    dat <- dat[, .SD, .SDcols=desired_cols]

    o <- foreach(chr=1:21, .combine='rbind') %do% {
        dat.sub <- dat[CHR==chr]
        dat.sub[, cM := pos_to_cM[[as.character(chr)]](BP)]
        dat.sub[, nextSNP := shift(BP, n=-1L), by=CHR]
        dat.sub[, nextSNP_cM := pos_to_cM[[as.character(chr)]](nextSNP)]
        dat.sub[, cM_between_SNPs := nextSNP_cM - cM]
        return(dat.sub)
    }
    return(o)
}

get_clusters <- function(DT, cM_threshold) {
    o <- copy(DT)
    o[, tf:= ifelse(cM_between_SNPs <= cM_threshold, TRUE, FALSE)]    # TRUE if next SNP is within cM_THRESHOLD, otherwise FALSE
    o[ is.na(nextSNP), tf := FALSE]

    # vectors of what cluster start and stop ROWS from o should be
    starts <- unique(c(1, which(o$tf == FALSE)+1))
    stops <- unique(c(which(o$tf == FALSE), nrow(o)))

    o2 <- foreach(start=starts, stop=stops, grp=1:length(starts), .combine='rbind') %do% {
        o.sub <- o[start:stop]
        o.sub[, 'cluster' := grp]
        return(o.sub[])
    }

    o2[, tf := NULL]

    o3 <- o2[, list('bp_start'=min(BP),
              'bp_stop'=max(BP),
              'cM_start'=min(cM),
              'cM_stop'=max(cM),
              'nSNPs_in_cluster'=.N), by=list(CHR,cluster)]
    
    o3[, 'cM_threshold' := cM_threshold]
    return(o3)
}

extend_clusters <- function(DT, cM_to_extend) {
    o <- foreach(chr=1:21, .combine='rbind') %do% {
        dat.sub <- copy(DT[CHR==chr])
        dat.sub[, extended_bp_start := cM_to_pos[[as.character(chr)]](cM_start - cM_to_extend)]
        dat.sub[, extended_bp_start := floor(extended_bp_start)]    # down to nearest integer
        dat.sub[is.na(extended_bp_start), extended_bp_start := 1]   # use start of chr if NA
        dat.sub[, extended_bp_stop := cM_to_pos[[as.character(chr)]](cM_start + cM_to_extend)]
        dat.sub[, extended_bp_stop := ceiling(extended_bp_stop)]    # up to nearest integer
        dat.sub[is.na(extended_bp_stop), extended_bp_stop := 3e8]   # use sufficiently large value to include end of chr if NA
        return(dat.sub)
    }
    return(o)
}


#### RUN ###########################################################################################

# Import recombination rates
bed <- import_recombination_bed(RECOMBINATION_BED_FILE)

# Build cM <-> BP functions
pos_to_cM <- build_pos_to_cM(bed)
cM_to_pos <- build_cM_to_pos(bed)

# Import gwas data
signif_snps <- get_gwas_hits(GWAS_FILENAME)
GWAS_NAME <- get_gwas_filestem(GWAS_FILENAME)
OUT_STEM=paste0(OUTDIR,'/',GWAS_NAME)

# Test a range of cM values
clusters <- foreach(cM=seq(0.01, 1, 0.01), .combine='rbind') %do% {
    get_clusters(signif_snps, cM)
}
clusters[, 'GWAS' := GWAS_NAME]

setkey(clusters, GWAS, cM_threshold, CHR, cluster)
desired_col_order <- c('GWAS', 'cM_threshold', 'cluster', 'CHR', 'bp_start', 'bp_stop',
                    'cM_start', 'cM_stop','nSNPs_in_cluster')

clusters <- clusters[, .SD, .SDcols=desired_col_order]
cluster_sizes <- clusters[, list('Clusters'=.N), by=list(cM_threshold)]
setkey(cluster_sizes, cM_threshold)

clusters_max <- max(cluster_sizes$Clusters)
clusters_min <-  min(cluster_sizes$Clusters)
cluster_N_range <-  clusters_max - clusters_min


# c_min <- cluster_sizes[which.max(cM_threshold)]$Clusters
# cluster_sizes[, ht := Clusters - c_min]
# ht_sum <- sum(cluster_sizes$ht)

cluster_sizes[, cumulative_frac := (clusters_max - Clusters)/cluster_N_range]

selected_cM_threshold <- cluster_sizes[cumulative_frac >= 0.9][1,]$cM_threshold


fwrite(clusters, file=paste0(OUT_STEM,'.clusters_all.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

g <- ggplot(cluster_sizes , aes(x=cM_threshold, y=Clusters)) +
    geom_point(shape=21) +
    theme_few(8) +
    scale_x_continuous(breaks=seq(0,1,0.2), labels=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(x='cM clustering threshold',
        y='Count',
        title='Clusters across all chromosomes') +
    ylim(0,NA)

ggsave(g, file=paste0(OUT_STEM, '.clusters.png'), width=10, height=10, units='cm')

g2 <- ggplot(cluster_sizes , aes(x=cM_threshold, y=cumulative_frac)) +
    geom_hline(yintercept=0.9, alpha=0.5, linetype='dashed') +
    geom_vline(xintercept=selected_cM_threshold, alpha=0.5, linetype='dashed', color='red') +
    geom_point(shape=21) +
    theme_few(8) +
    scale_x_continuous(breaks=sort(unique(c(1,selected_cM_threshold)))) +
    labs(x='cM clustering threshold',
        y='Effective Reduction in Cluster Count') +
    scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,NA))

ggsave(g2, file=paste0(OUT_STEM, '.clusters_frac.png'), width=10, height=10, units='cm')


clusters <- extend_clusters(clusters[cM_threshold==selected_cM_threshold], selected_cM_threshold/2)

fwrite(clusters, file=paste0(OUT_STEM, '.clusters_chosen.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

quit(status=0)


