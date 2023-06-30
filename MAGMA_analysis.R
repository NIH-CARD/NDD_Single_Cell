# Script to run MAGMA
# Load in libraries
library(MAGMA.Celltyping) 
library(dplyr)
library(EWCE)
library(knitr)
library(ggplot2)

# paths to formatted_sumstats
ad_sum <- "~/Documents/single_cell/AD_Bellenguez.formatted.tsv" # weird path set up needed? don't include MAGMA_Files in path only everything before and the original sumstats name?
lbd_sum <- "~/Documents/single_cell/LBD_Chia.formatted.tsv"
pd_sum <- "~/Documents/single_cell/PD_Nalls.formatted.tsv"
als_sum <- "~/Documents/single_cell/ALS_vanRheenen.formatted.tsv"
psp_sum <-"~/Documents/single_cell/PSP_Hoglinger.formatted.tsv"
ftld_sum <- "~/Documents/single_cell/FTLD_Pottier.formatted.tsv"

# path to ctd file
ctd_path <- "~/Documents/single_cell/ctd_adult_cluster_genes_mod.rda" 

# load in CTD
ctd <- load_rdata(ctd_path)

# Run the main cell type association analysis

# AD
ctAssocsLinear_AD <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = ad_sum, 
  analysis_name = "Linear",
  ctd_species = "human",
  force_new = TRUE
  )

FigsLinear_AD <- plot_celltype_associations(
  ctAssocs = ctAssocsLinear_AD,
  ctd = ctd)

ctAssocsTop_AD <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = ad_sum,
  analysis_name = "Top10",
  ctd_species = "human",
  EnrichmentMode = "Top 10%"#,
  #force_new = TRUE
  )

FigsTopDecile_AD <- plot_celltype_associations(
  ctAssocs = ctAssocsTop_AD,
  ctd = ctd)

ctAssocMerged_AD <- merge_magma_results(
  ctAssoc1 = ctAssocsLinear_AD,
  ctAssoc2 = ctAssocsTop_AD)

FigsMerged_AD <- plot_celltype_associations(
  ctAssocs = ctAssocMerged_AD,
  ctd = ctd)

# level 1
results_ad <- ctAssocMerged_AD[[1]]$results
tile_ad <- magma_tileplot(ctd=ctd, height = 10 , width = 10, results=results_ad, plotDendro=NULL, bind_plots=NULL, output_path = '~/Documents/single_cell/tileplots/')
ggsave(tile_ad, file='AD_tile_l1_dendro.png', width=20, height=20, units='cm', dpi=600)

# level 2
results_ad_l2 <- ctAssocMerged_AD[[2]]$results
tile_ad2 <- magma_tileplot(ctd=ctd, results=results_ad_l2, annotLevel = 2)
ggsave(tile_ad2, file='AD_tile_l2.png', width=20, height=20, units='cm', dpi=600)

# ALS
ctAssocsLinear_als_gwas <- calculate_celltype_associations(
  ctd = ctd,
  analysis_name = "Linear",
  gwas_sumstats_path = als_sum, 
  ctd_species = "human"#,
  #force_new = TRUE
  )

FigsLinear_als_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocsLinear_als_gwas,
  ctd = ctd)

ctAssocsTop_als_gwas <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = als_sum, 
  ctd_species = "human",
  EnrichmentMode = "Top 10%"#,
  #force_new = TRUE
  )

FigsTopDecile_als_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocsTop_als_gwas,
  ctd = ctd)

ctAssocMerged_als_gwas <- merge_magma_results(
  ctAssoc1 = ctAssocsLinear_als_gwas,
  ctAssoc2 = ctAssocsTop_als_gwas)

FigsMerged_als_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocMerged_als_gwas,
  ctd = ctd)

# level 1
results_als_l1 <- ctAssocMerged_als_gwas[[1]]$results
tile_als1 <- magma_tileplot(ctd=ctd, results=results_als_l1)

# level 2
results_als_l2 <- ctAssocMerged_als_gwas[[2]]$results
tile_als2 <- magma_tileplot(ctd=ctd, results=results_als_l2, annotLevel = 2)

# FTLD
ctAssocsLinear_ftld_gwas <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = ftld_sum,
  analysis_name = "Linear",
  ctd_species = "human"#,
  #force_new = TRUE
  )

FigsLinear_ftld_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocsLinear_ftld_gwas,
  ctd = ctd)

ctAssocsTop_ftld_gwas <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = ftld_sum, 
  ctd_species = "human",
  EnrichmentMode = "Top 10%"#,
  #force_new = TRUE
  )

FigsTopDecile_ftld_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocsTop_ftld_gwas,
  ctd = ctd)

ctAssocMerged_ftld_gwas <- merge_magma_results(
  ctAssoc1 = ctAssocsLinear_ftld_gwas,
  ctAssoc2 = ctAssocsTop_ftld_gwas)

FigsMerged_ftld_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocMerged_ftld_gwas,
  ctd = ctd)

# level 1
results_ftld <- ctAssocMerged_ftld_gwas[[1]]$results
tile_ftld1 <- magma_tileplot(ctd=ctd, results=results_ftld)

# level 2
results_ftld_l2 <- ctAssocMerged_ftld_gwas[[2]]$results
tile_ftld2 <- magma_tileplot(ctd=ctd, results=results_ftld_l2, annotLevel = 2, qvalue_thresh = 0.05)

## LBD
ctAssocsLinear_LBD <- calculate_celltype_associations(
  ctd = ctd,
  analysis_name = "Linear",
  gwas_sumstats_path = lbd_sum, 
  ctd_species = "human",
  force_new = TRUE)

FigsLinear_LBD <- plot_celltype_associations(
  ctAssocs = ctAssocsLinear_LBD,
  ctd = ctd)

ctAssocsTop_LBD <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = lbd_sum, 
  ctd_species = "human",
  EnrichmentMode = "Top 10%",
  force_new = TRUE)

FigsTopDecile_LBD <- plot_celltype_associations(
  ctAssocs = ctAssocsTop_LBD,
  ctd = ctd)


ctAssocMerged_LBD <- merge_magma_results(
  ctAssoc1 = ctAssocsLinear_LBD,
  ctAssoc2 = ctAssocsTop_LBD)


FigsMerged_LBD <- plot_celltype_associations(
  ctAssocs = ctAssocMerged_LBD,
  ctd = ctd)

# level 1
results_lbd <- ctAssocMerged_LBD[[1]]$results
tile_lbd1 <- magma_tileplot(ctd=ctd, results=results_lbd)
ggsave(tile_lbd1, file='LBD_tile_l1.png', width=20, height=20, units='cm', dpi=600)

# level 2
results_lbd_l2 <- ctAssocMerged_LBD[[2]]$results
tile_lbd2 <- magma_tileplot(ctd=ctd, results=results_lbd_l2, annotLevel = 2)
ggsave(tile_lbd2, file='LBD_tile_l2.png', width=20, height=20, units='cm', dpi=600)

## PD
ctAssocsLinear_PD <- calculate_celltype_associations(
  ctd = ctd,
  analysis_name = "Linear",
  gwas_sumstats_path = pd_sum, 
  ctd_species = "human",
  #force_new = TRUE
  )

FigsLinear_PD<- plot_celltype_associations(
  ctAssocs = ctAssocsLinear_PD,
  ctd = ctd)

ctAssocsTop_PD <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = pd_sum, 
  ctd_species = "human",
  EnrichmentMode = "Top 10%",
  #force_new = TRUE
  )

FigsTopDecile_PD <- plot_celltype_associations(
  ctAssocs = ctAssocsTop_PD,
  ctd = ctd)

ctAssocMerged_PD <- merge_magma_results(
  ctAssoc1 = ctAssocsLinear_PD,
  ctAssoc2 = ctAssocsTop_PD)

FigsMerged_PD <- plot_celltype_associations(
  ctAssocs = ctAssocMerged_PD,
  ctd = ctd)

# level 1
results_pd <- ctAssocMerged_PD[[1]]$results
tile_pd1 <- magma_tileplot(ctd=ctd, results=results_pd)
ggsave(tile_pd1, file='PD_tile_l1.png', width=20, height=20, units='cm', dpi=600)

# level 2
results_pd_l2 <- ctAssocMerged_PD[[2]]$results
tile_pd2 <- magma_tileplot(ctd=ctd, results=results_pd_l2, annotLevel = 2)
ggsave(tile_pd2, file='PD_tile_l2.png', width=20, height=20, units='cm', dpi=600)


# PSP
ctAssocsLinear_psp_gwas <- calculate_celltype_associations(
  ctd = ctd,
  analysis_name = "Linear",
  gwas_sumstats_path = psp_sum, 
  ctd_species = "human",
  force_new = TRUE)

FigsLinear_psp_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocsLinear_psp_gwas,
  ctd = ctd)

ctAssocsTop_psp_gwas <- calculate_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = psp_sum, 
  ctd_species = "human",
  EnrichmentMode = "Top 10%",
  force_new = TRUE)

FigsTopDecile_psp_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocsTop_psp_gwas,
  ctd = ctd)

ctAssocMerged_psp_gwas <- merge_magma_results(
  ctAssoc1 = ctAssocsLinear_psp_gwas,
  ctAssoc2 = ctAssocsTop_psp_gwas)

FigsMerged_psp_gwas <- plot_celltype_associations(
  ctAssocs = ctAssocMerged_psp_gwas,
  ctd = ctd)


# Conditional Analysis
ctCondAssocs_AD <- MAGMA.Celltyping::calculate_conditional_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = ad_sum, 
  analysis_name = "Conditional",
  controlTopNcells= 1)

ctCondAssocs_LBD <- calculate_conditional_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = lbd_sum, 
  analysis_name = "Conditional",
  controlTopNcells= 1,
  controlledAnnotLevel = 2,
  force_new = TRUE)

ctCondAssocs_PD <- calculate_conditional_celltype_associations(
  ctd = ctd,
  gwas_sumstats_path = pd_sum, 
  analysis_name = "Conditional",
  controlledCTs = c("none"),
  controlledAnnotLevel = 2,
  force_new = TRUE)

cond_fig_pd <- plot_celltype_associations(
  ctAssocs = ctCondAssocs_PD,
  ctd = ctd)
