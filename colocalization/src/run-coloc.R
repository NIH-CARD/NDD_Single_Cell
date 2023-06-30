#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(doMC)
library(parallel)
registerDoMC(cores=4)
library(ggplot2)
library(ggthemes)
library(viridis)
library(coloc)

args <- commandArgs(trailingOnly=TRUE) #

TISSUE_eQTL_FILENAME <- args[1]     # TISSUE_eQTL_FILENAME <- 'data/eQTL/Cortex-EUR.signif.tsv.gz'
NDD_GWAS_FILENAME <- args[2]        # NDD_GWAS_FILENAME <- 'data/NDD/PD_Nalls.signif.tsv'
NDD_CLUSTERS_FILENAME <- args[3]    # NDD_CLUSTERS_FILENAME <- 'data/NDD/clusters/PD_Nalls.clusters_chosen.tsv'

TISSUE <- strsplit(basename(TISSUE_eQTL_FILENAME), split='\\.signif')[[1]][1]
NDD <- strsplit(basename(NDD_GWAS_FILENAME), split='\\.signif')[[1]][1]

# Import and format eQTL GWAS table
eQTL_GWAS <- fread(TISSUE_eQTL_FILENAME)
setnames(eQTL_GWAS, toupper(colnames(eQTL_GWAS)))
setkey(eQTL_GWAS, CHR, BP)


# Import and format NDD Cluster coordinate table
COLOC_CLUSTERS <- fread(NDD_CLUSTERS_FILENAME)
setkey(COLOC_CLUSTERS, GWAS, CHR, cluster)

# Subset eQTLs to only regions of interest
eQTL_GWAS <- foreach(CHROM=COLOC_CLUSTERS$CHR,
        start=COLOC_CLUSTERS$extended_bp_start,
        stop=COLOC_CLUSTERS$extended_bp_stop,
        .combine='rbind', .errorhandling='remove') %do% {
            eQTL_GWAS[.(CHROM)][BP %between% c(start, stop)]
}
setkey(eQTL_GWAS, CHR, BP)



# Import and format NDD GWAS table
NDD_GWAS <- fread(NDD_GWAS_FILENAME)
setnames(NDD_GWAS, toupper(colnames(NDD_GWAS)))
setkey(NDD_GWAS, CHR, BP)

# Ensure eQTL alleles are polarized the same as NDD
NDD_GWAS[, 'chr_snp' := paste0(CHR, '_', BP, '_', A1, '_', A2)]
eQTL_GWAS[, 'chr_snp_rev' := paste0(CHR, '_', BP, '_', A2, '_', A1)]

eQTL_GWAS[chr_snp_rev %in% NDD_GWAS$chr_snp, 'rev' := TRUE]
eQTL_GWAS[rev==TRUE, B := -1*B]          # convert beta to other allele
eQTL_GWAS[rev==TRUE, FREQ := 1-FREQ]     # convert freq to other allele
eQTL_GWAS[rev==TRUE, a1 := A2]                               # rename
eQTL_GWAS[rev==TRUE, a2 := A1]                               # rename
eQTL_GWAS[rev==TRUE, A1 := a1]              # rename
eQTL_GWAS[rev==TRUE, A2 := a2]              # rename
eQTL_GWAS[, a1 := NULL]
eQTL_GWAS[, a2 := NULL]
eQTL_GWAS[, rev := NULL]
eQTL_GWAS[, chr_snp_rev := NULL]

eQTL_GWAS[, snpid := tstrsplit(SNP, split=':')[3]]
eQTL_GWAS[, SNP := paste0(CHR, ':', BP, ':', snpid, ':', A1, '_', A2)]
NDD_GWAS[, SNP := paste0(CHR, ':', BP, ':', SNP, ':', A1, '_', A2)]
eQTL_GWAS[, snpid := NULL]
NDD_GWAS[, chr_snp := NULL]

eQTL_GWAS[, 'id' := paste0(CHR, '_', BP, '_', PROBE)]
eQTL_GWAS[, 'rank' := frank(P), by=list(CHR,BP)]

gc()


run_coloc <- function(eQTL.DT, diagnosis.DT) {
    d <- copy(diagnosis.DT)
    diagnosis <- list()
    diagnosis$position <- d$BP
    diagnosis$snp <- d$SNP
    diagnosis$beta <- d$BETA
    diagnosis$varbeta <- d$SE^2
    diagnosis$type <- 'cc'

    e <- copy(eQTL.DT)
    eQTL <- list()
    eQTL$position <- e$BP
    eQTL$snp <- e$SNP
    eQTL$beta <- e$B
    eQTL$varbeta <- e$SE^2
    eQTL$type <- 'quant'
    eQTL$sdY <- rep(1, length(eQTL$snp))

    res <- coloc.abf(diagnosis, eQTL)
    res <- res$results
    setDT(res)
    return(res)
}



iteration <- 0

# Init first run

iteration <- iteration + 1
dt.all <- foreach(CHROM=COLOC_CLUSTERS$CHR,
        start=COLOC_CLUSTERS$extended_bp_start,
        stop=COLOC_CLUSTERS$extended_bp_stop,
        .combine='rbind', .errorhandling='remove') %do% {
            ndd_gwas <- NDD_GWAS[.(CHROM)][BP %between% c(start, stop)]
            eQTL_gwas <- eQTL_GWAS[.(CHROM)][rank==1 & BP %between% c(start, stop)]
            coloc.dt <- run_coloc(eQTL_gwas, ndd_gwas)
            coloc.dt[, 'iteration' := iteration]
            return(coloc.dt)
}

if(is.null(dt.all)) {quit(status=0)}

dt.all <- unique(dt.all)
dt.all <- dt.all[SNP.PP.H4 > 0.9]
dt.all[, c('CHR','BP') := tstrsplit(snp, split=':')[1:2]]
dt.all[, CHR := as.numeric(CHR)]
dt.all[, BP := as.numeric(BP)]
setkey(dt.all, CHR, BP)
all_coloc_results <- dt.all



d1 <- unique(dt.all[, .SD, .SDcols=c('CHR','BP')])
d1 <- unique(merge(d1, eQTL_GWAS[rank==1]))
d1[, 'iteration' := iteration]
all_d1 <- d1

eQTL_GWAS <- eQTL_GWAS[! id %in% d1$id]
eQTL_GWAS[,'rank' := frank(P), by=list(CHR, BP)]    # Reset eQTL probe rank





while(TRUE) {
# Iterate through 6 times, at which point no more signif
    iteration <- iteration + 1
    dt.all <- foreach(CHROM=COLOC_CLUSTERS$CHR,
            start=COLOC_CLUSTERS$extended_bp_start,
            stop=COLOC_CLUSTERS$extended_bp_stop,
            .combine='rbind', .errorhandling='remove') %do% {
                diagnosis_gwas <- NDD_GWAS[.(CHROM)][BP %between% c(start, stop)]
                eQTL_gwas <- eQTL_GWAS[.(CHROM)][rank==1 & BP %between% c(start, stop)]
                coloc.dt <- run_coloc(eQTL_gwas, diagnosis_gwas)
                coloc.dt[, 'iteration' := iteration]
                return(coloc.dt)
    }
    if(is.null(dt.all)) {break}
    dt.all <- unique(dt.all)
    dt.all <- dt.all[SNP.PP.H4 > 0.9]
    if(nrow(dt.all)==0) {break}
    dt.all[, c('CHR','BP') := tstrsplit(snp, split=':')[1:2]]
    dt.all[, CHR := as.numeric(CHR)]
    dt.all[, BP := as.numeric(BP)]
    setkey(dt.all, CHR, BP)
    all_coloc_results <- rbindlist(list(all_coloc_results, dt.all))



    d1 <- unique(dt.all[, .SD, .SDcols=c('CHR','BP')])
    d1 <- unique(merge(d1, eQTL_GWAS[rank==1]))
    d1[, 'iteration' := iteration]
    all_d1 <- rbindlist(list(all_d1, d1))
    eQTL_GWAS <- eQTL_GWAS[! id %in% d1$id]
    eQTL_GWAS[,'rank' := frank(P), by=list(CHR, BP)]    # Reset eQTL probe rank
}



setkey(all_coloc_results, iteration, CHR, BP)
setkey(all_d1, iteration, CHR, BP)

dt.out <- merge(all_coloc_results, all_d1)
dt.out[, 'TISSUE' := TISSUE]
dt.out[, 'NDD' := NDD]
dt.out[, snp := NULL]

fwrite(dt.out, file=paste0('data/coloc/', NDD, '-', TISSUE, '.tsv'), quote=F, row.names=F, col.names=T, sep='\t')

quit(status=0)

















###

d2 <- unique(eQTL_GWAS[rank==1, .SD, .SDcols=c(desired_cols,'PROBE')])
d2[, id := paste0(CHR,'-',BP,'-',TISSUE)]

d2[! id %in% d1$id]
d2[ id %in% d1$id]

fsetdiff(d1, d2)

fwrite(dt.all, file='data/coloc_ALL.tsv', row.names=F, col.names=T, sep='\t', quote=FALSE)

dt.formatted <- dt.all[, .SD, .SDcols=c('diagnosis','TISSUE','rank','CHR','BP','SNP.PP.H4')]
setnames(dt.formatted, 'SNP.PP.H4', 'PP')

dt.signif <- dt.formatted[PP > 0.90]

setkey(dt.signif, TISSUE, CHR, BP, rank)

eQTL_tomerge <- eQTL_GWAS[, .SD, .SDcols=c('TISSUE','CHR','BP','PROBE','GENE','B','rank')]
setkey(eQTL_tomerge, TISSUE, CHR, BP, rank)


dt.merged <- merge(dt.signif, eQTL_tomerge)
setkey(dt.merged, 'CHR', 'BP', 'rank')
dt.merged[, 'rank' := NULL]

setnames(dt.merged, toupper(colnames(dt.merged)))
dt.merged[, GENE := NULL]

dt.merged <- dt.merged[, .SD, .SDcols=c('CHR','BP','DIAGNOSIS','TISSUE','PROBE','B','PP')]


fwrite(dt.merged, file='coloc_PP.tsv', quote=F, row.names=F, col.names=T, sep='\t')

quit(status=0)

dt.wide <- dcast(dt.merged, CHR+BP+diagnosis+TISSUE~probe_rank, value.var='PROBE')

fwrite(dt.wide, file='coloc_top3_probes_wide.tsv', quote=F, row.names=F, col.names=T, sep='\t')






# Take the most significant probe (or gene, same result) per SNP
# Assumption that SNP influences only one gene
# eQTLs <- eQTLs[eQTLs[, .I[which.min(P)], by=list(CHR,BP)]$V1]

# Sort on keys
setkey(eQTLs, CHR, BP)
# function the same as:
# eQTLs <-  eQTLs[order(CHR,BP)]


eQTLs[, 'TISSUE' := TISSUE]

diagnosis_clusters <- copy(all_clusters[GWAS==diagnosis
























DISEASE_GWAS_COLS <- c('SNP','CHR','BP','A1','A2','BETA','FRQ','SE','P')
eQT_GWAS_COLS <- 

# Read in .bed file containing recombination rates
# Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
bed <- fread(file = RECOMBINATION_BED_FILE,
              header = TRUE,
              showProgress = FALSE,
              col.names = c("CHR","BP","c","cM")
)
bed[, c := NULL]

setkey(bed, CHR, BP)
bed_extension <- bed[, list('BP'=max(BP)*1.1, 'cM'=max(cM)), by=CHR]
bed <- rbindlist(list(bed, bed_extension))
setkey(bed, CHR, BP)


# Generate functions (cM ~ BP) to translate base pair position to units of cM
cM_to_POS <- new.env()
POS_to_cM <- new.env()
for(chr in unique(bed[,CHR])) {
    cM_to_POS[[as.character(chr)]] <- approxfun(y=c(0, bed[CHR==chr][,BP]), x=c(0,bed[CHR==chr][,cM]), ties='ordered')
    POS_to_cM[[as.character(chr)]] <- approxfun(x=c(0, bed[CHR==chr][,BP]), y=c(0,bed[CHR==chr][,cM]), ties='ordered')
}





# get_linkage_window <- function(DT.disease, DT.eQTL, DT.linkage, chr, pos) {
#     setnames(DT.linkage, c('CHR','POS','RATE','cM'))
#     original_cM <- DT.linkage[CHR==chr & POS==pos][['position']]
#     cM_window <- DT.linkage[CHR==chr & POS %between% c((pos-1), (pos+1))][['POS']]
#     min_cM <- min(cM_window)
#     max_cM <- max(cM_window)
#     return(cM_window)
#     #return(c(min_cm, max_cm))
# }

# get_linkage_window(linkage, linkage, linkage, 1, 61551035)


For colocalization analyses, we included data from eQTL and disease GWAS within one map unit (centimorgan, cM) of a 

map units


# Read in .bed file containing recombination rates
# Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
bed <- fread(file = RECOMBINATION_BED_FILE,
              header = TRUE,
              showProgress = FALSE,
              col.names = c("chr","pos","c","cM")
)

# bed <- bed[chr==CHROM]
setkey(bed, chr, pos)

# take latest value of cM
# bed[ !duplicated(bed[, c("cM")], fromLast=T),][order(chr, pos, cM)]


# Generate functions (cM ~ BP) to translate base pair position to units of cM
cM_to_POS <- new.env()
POS_to_cM <- new.env()

# Create recombination function based on .bed file
for(chr.i in unique(bed[,chr])) {
    cM_to_POS[[as.character(chr.i)]] <- approxfun(y=c(0, bed[chr==chr.i][,pos]), x=c(0,bed[chr==chr.i][,cM]), ties='ordered')
    POS_to_cM[[as.character(chr.i)]] <- approxfun(x=c(0, bed[chr==chr.i][,pos]), y=c(0,bed[chr==chr.i][,cM]), ties='ordered')
}


cM_to_POS[['1']](47.29347)
POS_to_cM[['1']](24892280)


Predicted positional cM by linearly interpolating map units from Broad Institute Hg38 linkage map. 

bed[chr==1, pos]
dt <- data.table(actual_cM = bed[chr==1,cM],
           predicted_cM =recombination_function[['1']](bed[chr==1, pos]))

ggplot(data=dt, aes(x=actual_cM, y=predicted_cM)) + geom_point() 

ggplot(bed, 
































get_eQTL <- function(datadir, tissue, chr) {
    eQTL_fn <- paste0(datadir, '/', tissue, '_', chr, '.tsv.gz')
    fread(eQTL_fn)
}

tissue <- 'Cortex-EUR'
chr <- 21

dat <- get_eQTL(eQTL_datadir, tissue, chr)




# read in formatted summary statistics - gr38

# AD_fn <- 'data/AD_Bellenguez.formatted.tsv.gz'
# LBD_fn <- 'data/LBD_Chia.formatted.tsv.gz'
# PD_fn <- 'data/PD_Nalls.formatted.tsv.gz'

# basalganglia_fn <- 'data/Basalganglia-EUR-out.cis.summary.txt'
# cerebellum_fn <- 'data/Cerebellum-out.cis.summary.txt'
# cortex_fn <- 'data/Cortex-EUR-out.cis.summary.txt'
# hippocampus_fn <- 'data/Hippocampus-out.cis.summary.txt'
# spinalcord_fn <- 'data/Spinalcord-out.cis.summary.txt'

diagnoses <- list(c('AD','data/AD_Bellenguez.formatted.tsv.gz'),
                    c('LBD','data/LBD_Chia.formatted.tsv.gz'),
                    c('PD','data/PD_Nalls.formatted.tsv.gz')
)

celltypes <- list(c('basal_ganglia','data/Basalganglia-EUR-out.cis.summary.txt'),
                c('cerebellum','data/Cerebellum-out.cis.summary.txt'),
                c('cortex','data/Cortex-EUR-out.cis.summary.txt'),
                c('hippocampus','data/Hippocampus-out.cis.summary.txt'),
                c('spinal_cord','data/Spinalcord-out.cis.summary.txt')
)

# coloc can specify sdY=1 if expression was normalized to standard normal distribution

compiled_table_fn <- 'merged_hits.tsv'

DT[, varbeta = se^2 ]

dat[, varbeta := se^2]

diagnoses <- c( 'AD'='data/AD_Bellenguez.formatted.tsv.gz',
                'LBD'='data/LBD_Chia.formatted.tsv.gz',
                'PD'='data/PD_Nalls.formatted.tsv.gz')

celltypes <- c( 'basal_ganglia'='data/Basalganglia-EUR-out.cis.summary.txt',
                'cerebellum'='data/Cerebellum-out.cis.summary.txt',
                'cortex'='data/Cortex-EUR-out.cis.summary.txt',
                'hippocampus'='data/Hippocampus-out.cis.summary.txt',
                'spinal_cord'='data/Spinalcord-out.cis.summary.txt')



get_diagnosis <- function(diagnosis_fn, diagnosis_name) {
    dat <- fread(diagnosis_fn)
    dat[, c('SNP','A1','A2') := NULL]
    setnames(dat, c('BETA','FRQ','SE','P'), c('beta','Freq', 'se', 'p'))
    dat[, 'diagnosis' := diagnosis_name]
    setkey(dat, CHR, BP)
    return(dat[])
}

get_eQTL <- function(celltype_fn, celltype_name) {
    dat <- fread(celltype_fn)
    dat <- dat[SNP_Chr != 'SNP_Chr']
    dat[,c('TopSNP','A1','A2','cis_Chr','cis_left','cis_right','Orientation','Probe_Chr','Probe', 'Probe_bp', 'Gene') := NULL]
    setnames(dat, 'SNP_Chr', 'CHR')
    dat <- dat[order(CHR,BP,p)]
    dat <- dat[which(!duplicated(dat[,c(1:2)]))]        # get non-duplicate
    dat[, CHR := as.numeric(CHR)]
    dat[, BP := as.numeric(BP)]
    dat[, p := as.numeric(p)]
    dat[, se := as.numeric(se)]
    dat[, b := as.numeric(b)]
    dat[, Freq := as.numeric(Freq)]
    setnames(dat, c('Freq','b','se','p'), c('Freq','beta','se','p'))
    setkey(dat, CHR, BP)
    dat[, 'celltype' := celltype_name]
    return(dat[])
}


run_coloc <- function(diagnosis.DT, eQTL.DT, chromosome) {
    d <- copy(diagnosis.DT)[CHR==chromosome]    # subset by chromosome
    diagnosis <- list()
    diagnosis$position <- d$BP
    diagnosis$snp <- paste0(chromosome, '_', d$BP)
    diagnosis$beta <- d$beta
    diagnosis$varbeta <- d$se^2
    diagnosis$type <- 'cc'


    e <- copy(eQTL.DT)[CHR==chromosome]
    eQTL <- list()
    eQTL$position <- e$BP
    eQTL$snp <- paste0(chromosome, '_', e$BP)
    eQTL$beta <- e$beta
    eQTL$varbeta <- e$se^2
    eQTL$type <- 'quant'
    eQTL$sdY <- rep(1, length(eQTL$snp))

    res <- coloc.abf(diagnosis, eQTL)
    res <- res$results
    setDT(res)
    return(res)
}

## ENSURE REF AND ALT ALLELE POLARIZED THE SAME!
o <- foreach(d_fn=diagnoses, d_name=names(diagnoses), .combine='rbind') %do% {
    diagnosis.DT <- get_diagnosis(d_fn, d_name)
    foreach(e_fn=celltypes, e_name=names(celltypes), .combine='rbind') %do% {
        eQTL.DT <- get_eQTL(e_fn, e_name)
        o1 <- foreach(chromosome=1:22, .combine='rbind') %do% {
            run_coloc(diagnosis.DT, eQTL.DT, chromosome)
        }
        o1[, 'celltype' := e_name]
        o1[, 'diagnosis' := d_name]
        return(o1[])
    }
}

saveRDS(o, file='coloc.RDS')

o.signif <- o[SNP.PP.H4 > 0.5, .SD, .SDcols=c('chr','position','diagnosis','celltype','SNP.PP.H4')]
o.signif <- o.signif[order(diagnosis, celltype, chr, position, SNP.PP.H4)]
fwrite(o.signif, file='coloc_50_PP.tsv', quote=F, row.names=F, col.names=T, sep='\t')


o[, chr := tstrsplit(snp, split='_')[1]]
o[, chr := as.numeric(chr)]


o2 <- order(o$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(o$SNP.PP.H4[o2])
w <- which(cs > 0.95)[1]
my.res$results[o,][1:w,]$snp

ggplot(o[SNP.PP.H4>0.95], aes(x=position, y=-1*log10(SNP.PP.H4), color=diagnosis, shape=celltype)) + geom_point() + facet_grid(chr~.) + theme_few()

subset(o,SNP.PP.H4>0.95)

check_dataset(eQTL)
check_dataset(diagnosis)

res <- coloc.abf(diagnosis, eQTL)

if (! file.exists(compiled_table_fn)) {
    dat <- foreach(i=diagnoses, .combine='rbind') %do% {
        diagnosis_type <- i[1]
        diagnosis_filename <- i[2]
        dt.diagnosis <- get_diagnosis(diagnosis_filename, diagnosis_type)
        foreach(j=celltypes, .combine='rbind') %do% {
            celltype_name <- j[1]
            celltype_filename <- j[2]
            dt.eQTL <- get_eQTL(celltype_filename, celltype_name)
            dt.merge <- merge(dt.diagnosis, dt.eQTL)
            return(dt.merge[])
        }
    }
    dat <- dat[eQTL_p < 0.05 & diagnosis_p < 0.05]
    #dat[, eQTL_hit := ifelse(eQTL_p < 0.05, TRUE, FALSE)]
    #dat[, diagnosis_hit := ifelse(diagnosis_p < 0.05, TRUE, FALSE)]

    # chrLengths <- data.table(CHR=1:22, bp=
    # c('248956422','242193529','198295559','190214555','181538259','170805979','159345973',
    # '145138636','138394717','133797422','135086622','133275309','114364328','107043718',
    # '101991189','90338345','83257441','80373285','58617616','64444167','46709983','50818468'))

    # chrLengths[, cumulative := shift(cumsum(bp), type='lag', fill=0)]
    # chrLengths[, bp := NULL]

    # dat  <- merge(dat, chrLengths, by.x='CHR', by.y='CHR')
    fwrite(dat, file='merged_hits.tsv', quote=F, row.names=F, col.names=T, sep='\t')
} else {
    dat <- fread(compiled_table_fn)
}


dat[, grp := paste0(celltype, '_', diagnosis)]


# Separated chromosomes

windowWidth <- 1e6
dat[, window := floor(BP / windowWidth), by=CHR]
dat.ag <- dat[, .N, by=list(diagnosis, celltype, window, CHR)]


allcombos <- foreach(celltype = sapply(celltypes, function(x) x[1]), .combine='rbind') %do% {
    foreach(diagnosis = sapply(diagnoses, function(x) x[1]), .combine='rbind') %do% {
        cbind(dat[, list('window'=0:max(window)), by=CHR][order(CHR, window)], 'celltype'=celltype, 'diagnosis'=diagnosis)
    }
}

setkey(dat.ag, celltype, diagnosis, CHR, window)
setkey(allcombos, celltype, diagnosis, CHR, window)

dat.merge <- merge(dat.ag, allcombos, all=TRUE)
dat.merge[is.na(N), N := 0]

for(chromosome in 1:22) {
    dat.sub <- dat.merge[CHR==chromosome]
    dat.sub[, N := factor(N)]
    outfile <- paste0(chromosome, '.png')

    g <- ggplot(data=dat.sub, aes(x=windowWidth*window, y=diagnosis, fill=N)) +
    geom_tile() +
    scale_fill_viridis(discrete=TRUE) +
    labs(x='500 kb window', y='') +
    theme_few() +
    facet_grid(celltype~CHR) +
    labs(title=paste0('chromosome ', chromosome), x=paste0('BP along chromosome, window width = ', windowWidth)) +
    theme(strip.text.x = element_blank(),
        strip.text.y = element_text(angle=0))

    ggsave(g, file=outfile, width=40, height=10, units='cm')

}

setkey(dat, CHR, BP)

# format GTF


split_keys <- function(x, all_keys) {
    split_line <- unlist(strsplit(x, split=';|; |;$'))
    out_cols <- c()
    for(key_value_pair in split_line) {
        split_pair <- unlist(strsplit(key_value_pair, split=' '))
        key_ <- gsub('"', '', split_pair[1])
        key_ <- gsub("'", '', key_)
        value_ <- gsub('"', '', split_pair[2])
        value_ <- gsub("'", '', value_)
        value_ <- unlist(value_)
        a <- data.frame('tmp'=value_)
        setnames(a, 'tmp', key_)
        out_cols <- c(out_cols, a)
    }
    out_cols <- as.data.table(rbind(out_cols))
    for(key_ in all_keys) {
        if(! key_ %in% colnames(out_cols)) {
            out_cols[, tmp := NA]
            setnames(out_cols, 'tmp', key_)
        }
    }
    setcolorder(out_cols, all_keys)
    return(out_cols)
}

#gtf_filename <- '/fdb/bcbio-nextgen/current/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.gtf'
gtf_filename <- '/fdb/ensembl/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf'
gtf <- fread(gtf_filename)
gtf <- gtf[V3=='gene']
setnames(gtf, colnames(gtf)[1], 'V1')
gtf <- gtf[V1 %in% 1:22]

gtf_0 <- copy(gtf)

gc()

all_keys <- unique(sapply(strsplit(unlist(strsplit(gtf$V9, split=';|; |;$')), split=' '), function(x) x[1]))

b <- foreach(i=gtf$V9, .combine='rbind') %do% {
    split_keys(i, all_keys)
}

gtf <- cbind(gtf, b)
gtf <- gtf[gene_biotype=='protein_coding']

if(length(unique(
            nrow(gtf),
            length(unlist(gtf$gene_name)),
            length(unlist(gtf$gene_id)),
            length(unlist(gtf$gene_version)),
            length(unlist(gtf$gene_source)),
            length(unlist(gtf$gene_biotype))
            )
) != 1) {
    cat('Row lengths for newly formed GTF are not equal!\n')
    quit(1)
}




gtf[, gene_name := unlist(gene_name)]
gtf[, gene_id := unlist(gene_id)]
gtf[, gene_version := unlist(gene_version)]
gtf[, gene_source := unlist(gene_source)]
gtf[, gene_biotype := unlist(gene_biotype)]
gtf[, V9 := NULL]


# subset to 21 autosomes, X, Y, mitochondrion
setnames(gtf, 'V1', 'CHR')
setnames(gtf, 'V4', 'start')
setnames(gtf, 'V5', 'stop')

gtf[, CHR := as.numeric(CHR)]

# Get SNPs inside genes
dat[, 'start' := BP]
dat[, 'stop' := BP]
dat[, idx := 1:.N]

setkey(gtf, CHR, start, stop)
overlaps <- foverlaps(dat, gtf, type="within")[!is.na(V3)]

noncoding_hits <- dat[! idx %in% unique(overlaps$idx)]
#noncoding_hits[, c('start','stop') := NULL]

# SO FAR SO GOOD

setkey(noncoding_hits, CHR, start)
setkey(gtf, CHR, start)



plus_search <- function(x, y) {
    # find the closest value in y where x comes before
    return(which.min(y<x))
}

minus_search <- function(x, y) {
    # find the closest value in y where x comes after
    return(max(which(y<x)))
}


o <- foreach(chr=1:22, .combine='rbind') %do% {
    dt.chr <- noncoding_hits[CHR==chr]
    dt.gtf.plus <- gtf[CHR==chr & V7 == '+']
    dt.gtf.minus <- gtf[CHR==chr & V7 == '-']
    plus_idx <- sapply(dt.chr$BP, function(x) plus_search(x, dt.gtf.plus$start))
    minus_idx <- sapply(dt.chr$BP, function(x) minus_search(x, dt.gtf.minus$stop))
    dt.gtf.plus[plus_idx]
    dt.gtf.minus[minus_idx]
    dt.chr[, plus_start := dt.gtf.plus[plus_idx]$start]
    dt.chr[, plus_gene_id := dt.gtf.plus[plus_idx]$gene_id]
    dt.chr[, plus_gene_name := dt.gtf.plus[plus_idx]$gene_name]
    
    dt.chr[, minus_start := dt.gtf.minus[minus_idx]$stop]
    dt.chr[, minus_gene_id := dt.gtf.minus[minus_idx]$gene_id]
    dt.chr[, minus_gene_name := dt.gtf.minus[minus_idx]$gene_name]

    dt.chr[, dist_to_plus_start := plus_start - BP]
    dt.chr[, dist_to_minus_start := BP - minus_start ]
    dt.chr[, 'closest_gene_id' := ifelse(dist_to_plus_start < dist_to_minus_start, plus_gene_id, minus_gene_id)]
    dt.chr[, 'closest_gene_name' := ifelse(dist_to_plus_start < dist_to_minus_start, plus_gene_name, minus_gene_name)]
}















want next MINUS position just below that
which.ma

want next PLUS position just after that
which.min(BP

exonic_hits[gtf, roll=Inf][!is.na(BP)]
exonic_hits[gtf, roll=-Inf][!is.na(BP)]

which.min = first TRUE
which.max = first FALSE

# Old overlaps
   # For faster computation, determine which ranges actually contain a SNP,
    # So that ranges without SNPs can be excluded downstream
        dat.sites <- copy(dat.vcf[,.(POS)])
        dat.sites[, nearest := POS]
        setkey(dat.sites, POS)
        setkey(segments.to.iter, start)
        segments.to.iter <- dat.sites[segments.to.iter, list(start, stop, nearest), roll=-Inf]
        segments.to.iter[, containsSNP := ifelse(nearest >= start & nearest <= stop, TRUE, FALSE)]


# edit ASE table to be foverlaps-able
dat.filtered[, "stop" := POS]
setnames(dat.filtered, "POS", "start")
setnames(dat.filtered, "id", "sampleID")
dat.filtered[, id := 1:.N]
setkey(dat.filtered, NULL)
# dat.wide.filtered[, "stop" := POS]
# setnames(dat.wide.filtered, "POS", "start")
# dat.wide.filtered[, id := 1:.N]



overlaps <- foverlaps(dat.filtered, exons, type="within")
overlaps[, id2 := rleid(id, GeneID)]
overlaps <- overlaps[!duplicated(id2) & !is.na(GeneID)]
overlaps[, c("start","stop","i.stop","id","id2") := NULL]
setnames(overlaps, "i.start", "POS")
setkey(overlaps, GeneID)
overlaps[, gn := rleid(GeneID)]






# Concatenated chromosomes
if(FALSE) {
dat[, cumulativeBP := BP + cumulative]
windowWidth <- 5e6
dat[, window := floor(cumulativeBP / windowWidth), by=CHR]
dat.ag <- dat[, .N, by=list(grp, window)]

allcombos <- CJ('grp'=unique(dat.ag$grp), 'window'=0:max(dat.ag$window))

dat.ag.all <- merge(dat.ag, allcombos, all=TRUE)
dat.ag.all[is.na(N), N:=0]

dat.ag.all[grp %like% "AD$", diagnosis := "AD"]
dat.ag.all[grp %like% "PD$", diagnosis := "PD"]
dat.ag.all[grp %like% "LBD$", diagnosis := "LBD"]

dat.ag.all[grp %like% "^basal_ganglia", celltype := "basal_ganglia"]
dat.ag.all[grp %like% "^spinal_cord", celltype := "spinal_cord"]
dat.ag.all[grp %like% "^cerebellum", celltype := "cerebellum"]
dat.ag.all[grp %like% "^cortex", celltype := "cortex"]
dat.ag.all[grp %like% "^hippocampus", celltype := "hippocampus"]
}

g <- ggplot(data=dat.ag.all, aes(x=window, y=diagnosis, fill=N)) +
geom_tile() +
scale_fill_viridis() +
labs(x='500 kb window', y='') +
theme_few() +
facet_grid(celltype~.)

g2 <- ggplot(data=dat, aes(x=window, y=diagnosis, fill=N)) +
geom_tile() +
scale_fill_viridis() +
labs(x='500 kb window', y='') +
theme_few() +
facet_grid(celltype~.)



ggsave(g, file='heatmap.png', width=60, height=8, units='cm')