#!/usr/bin/env Rscript

library(data.table)
library(foreach)

clusters.extended <- fread('data/clusters_extended.tsv')

DIAGNOSIS_GWAS_FILES <- list.files(path='/data/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/data/final_formatted_sumstats/',
                                    pattern='*.formatted.tsv',
                                    full.names=T)


get_gwas_filestem <- function(filename) {
    return(strsplit(basename(filename), split='\\.')[[1]][1])
}

diagnosis_gwas_SNPs <- foreach(filename = DIAGNOSIS_GWAS_FILES, .combine='rbind') %do% {
    diagnosis <- get_gwas_filestem(filename)
    dt <- fread(filename)  # already correct chromosome
    setnames(dt, toupper(colnames(dt)))
    desired_cols <- c('SNP','CHR','BP','A1','A2','BETA','FRQ','SE','P')
    dt <- dt[, .SD, .SDcols=desired_cols]

    foreach(CHROM=1:21, .combine='rbind') %do% {
        clusters <- clusters.extended[CHR==CHROM]
        o <- foreach(start=clusters$extended_bp_start, stop=clusters$extended_bp_stop, .combine='rbind') %do% {
            return(dt[CHR==CHROM & BP %between% c(start, stop)])      # Subset to window of interest
        }
        o <- unique(o)
        o[, 'DIAGNOSIS' := diagnosis]
        return(o)
    }
}

fwrite(diagnosis_gwas_SNPs, file='data/diagnosis_colocalization_SNPs.tsv', row.names=F, col.names=T, sep='\t', quote=FALSE)

