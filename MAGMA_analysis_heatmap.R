# Script to run MAGMA
# Load in libraries
library(MAGMA.Celltyping) 
library(dplyr)
library(EWCE)
library(knitr)

# paths to formatted_sumstats
ad_sum <- "~/Documents/single_cell/AD_Bellenguez.formatted.tsv" # weird path set up needed? don't include MAGMA_Files in path only everything before and the original sumstats name?
lbd_sum <- "~/Documents/single_cell/LBD_Chia.formatted.tsv"
pd_sum <- "~/Documents/single_cell/PD_Nalls.formatted.tsv"
als_sum <- "~/Documents/single_cell/ALS_vanRheenen.formatted.tsv"
psp_sum <-"~/Documents/single_cell/PSP_Hoglinger.formatted.tsv"
ftld_sum <- "~/Documents/single_cell/FTLD_Pottier.formatted.tsv"

# path to ctd file
ctd_path <- "~/Documents/single_cell/ctd_adult_cluster_genes.rda" 

# load in CTD
ctd <- load_rdata(ctd_path)

# Run the main cell type association analysis pipeline
## AD
ad <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = ad_sum,
  ctd = ctd,
  ctd_name = 'linnarson',
  ctd_species = "Human",
  run_linear = TRUE, 
  run_top10 = TRUE)
merged_AD <- MAGMA.Celltyping::merge_results(
  MAGMA_results = ad)
heat_AD_merge <- MAGMA.Celltyping::results_heatmap(
  merged_results = merged_AD , 
  title = "Alzheimer's Disease Bellenguez 2022 vs. Adult Human Brain Single Cells",
  fdr_thresh = 1)

tile_ad <- magma_tileplot(ctd,merged_AD)
## ALS
als <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = als_sum,
  ctd = ctd,
  ctd_name = 'linnarson',
  ctd_species = "Human",
  run_linear = TRUE, 
  run_top10 = TRUE)
merged_als <- MAGMA.Celltyping::merge_results(
  MAGMA_results = als)
heat_als_merge <- MAGMA.Celltyping::results_heatmap(
  merged_results = merged_als , 
  title = "Amyotrophic lateral sclerosis (ALS) van Rheenan 2021 vs. Adult Human Brain Single Cells",
  fdr_thresh = 1)

tile_als <- magma_tileplot(ctd,merged_als)
tile_als

## FTLD
ftld <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = ftld_sum,
  ctd = ctd,
  ctd_name = 'linnarson',
  ctd_species = "Human",
  run_linear = TRUE, 
  run_top10 = TRUE)
merged_ftld <- MAGMA.Celltyping::merge_results(
  MAGMA_results = ftld)
heat_ftld_merge <- MAGMA.Celltyping::results_heatmap(
  merged_results = merged_ftld , 
  title = "Frontotemporal lobar degeneration (FTLD) Pottier 2019 vs. Adult Human Brain Single Cells",
  fdr_thresh = 1)

## LBD
lbd <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = lbd_sum,
  ctd = ctd,
  ctd_name = 'linnarson',
  ctd_species = "Human",
  run_linear = TRUE, 
  run_top10 = TRUE)
merged_lbd <- MAGMA.Celltyping::merge_results(
  MAGMA_results = lbd)
heat_lbd_merge <- MAGMA.Celltyping::results_heatmap(
  merged_results = merged_lbd , 
  title = "Lewy body dementia (LBD) Chia 2021 vs. Adult Human Brain Single Cells",
  fdr_thresh = 1)

## PD
pd <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = pd_sum,
  ctd = ctd,
  ctd_name = 'linnarson',
  ctd_species = "Human",
  run_linear = TRUE, 
  run_top10 = TRUE)
merged_pd <- MAGMA.Celltyping::merge_results(
  MAGMA_results = pd)
heat_pd_merge <- MAGMA.Celltyping::results_heatmap(
  merged_results = merged_pd, 
  title = "Parkinson's disease (PD) Nalls 2019 vs. Adult Human Brain Single Cells",
  fdr_thresh = 1)

## PSP
psp <- MAGMA.Celltyping::celltype_associations_pipeline(
  magma_dirs = psp_sum,
  ctd = ctd,
  ctd_name = 'linnarson',
  ctd_species = "Human",
  run_linear = TRUE, 
  run_top10 = TRUE)
merged_psp <- MAGMA.Celltyping::merge_results(
  MAGMA_results = psp)
heat_psp_merge <- MAGMA.Celltyping::results_heatmap(
  merged_results = merged_psp , 
  title = "Progressive supranuclear palsy (PSP) Hoglinger 2011 vs. Adult Human Brain Single Cells",
  fdr_thresh = 1)