#!/usr/bin/env Rscript

library(data.table)

disease_gwas_files <- list.files(path='/data/CARD_AA/projects/2022_10_CA_singlecell_humanbrain/data/final_formatted_sumstats', full.names=T)
signif_threshold <- 5e-8

for(fn in disease_gwas_files) {
    dat <- fread(fn)
    base_filename <- basename(fn)
    gwas_name <- strsplit(base_filename, split='\\.formatted')[[1]][1]
    out_filename <- paste0(gwas_name, '.signif.tsv')

    dat <- dat[P <= signif_threshold]
    fwrite(dat, file=paste0('data/', out_filename), quote=F, row.names=F, col.names=T, sep='\t')
}