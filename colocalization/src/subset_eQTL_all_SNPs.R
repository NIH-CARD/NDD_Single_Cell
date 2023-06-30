#!/usr/bin/env Rscript

library(data.table)
library(foreach)

clusters.extended <- fread('data/clusters_extended.tsv')


for(CHROM in 1:21) {
    clusters <- clusters.extended[CHR==CHROM]
    filenames <- list.files(path='/data/CARD/projects/singlecell_humanbrain/coloc_20230322',
                            pattern=paste0('_', CHROM, '.tsv.gz'),
                            full.names=T)

    chr_eQTL_snps <- foreach(filename=filenames, .combine='rbind') %do% {
        tissue <- strsplit(basename(filename), split='_')[[1]][1]
        dt <- fread(filename)  # already correct chromosome
        setnames(dt, toupper(colnames(dt)))
        o <- foreach(start=clusters$extended_bp_start, stop=clusters$extended_bp_stop, .combine='rbind') %do% {
            dt <- dt[CHR==CHROM & BP %between% c(start, stop)]      # Subset to window of interest
            #dt <- dt[dt[, .I[which.min(P)], by=list(CHR,BP)]$V1]    # Take most significant probe per SNP
            return(dt)
        }
        o <- unique(o)
        o[, 'TISSUE' := tissue]
        return(o)
    }
    fwrite(chr_eQTL_snps, file=paste0('data/', CHROM, '_eQTL_all_SNPs.tsv'), row.names=F, col.names=T, sep='\t', quote=FALSE)
}

# Combining tables

dat <- foreach(file=list.files(path='data/', pattern='*_eQTL_all_SNPs.tsv', full.names=T), .combine='rbind') %do% {
    fread(file)
}

fwrite(dat, file='data/eQTL_colocalization_all_SNPs.tsv', quote=F, row.names=F, col.names=T, sep='\t')