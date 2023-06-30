# necessary installs
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38", force = TRUE) # download GR38 data
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37", force = TRUE) # download GR37 data
BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

# library load in
library(MAGMA.Celltyping) 
library(dplyr)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(MungeSumstats)

# Prepare Data
## GWAS

# paths to sum stats
AD_gwas = "/Users/alvaradocx/Documents/single_cell/AD_Bellenguez_sumstats.txt"
PD_gwas = "/Users/alvaradocx/Documents/single_cell/PD_nalls_sumstats.txt"
PD_meta_gwas = "/Users/alvaradocx/Documents/single_cell/PD_META_sumstats.txt"
LBD_gwas = "/Users/alvaradocx/Documents/single_cell/LBD_Chia_sumstats.txt"
PSP_gwas = "/Users/alvaradocx/Documents/single_cell/PSP_sumstats37.txt"
FTLD_gwas = "/Users/alvaradocx/Documents/single_cell/FTLD_GWAS37.csv"
FTLD2_gwas = "/Users/alvaradocx/Documents/single_cell/FTLD_gwas38_test.txt"
ALS_gwas = "/Users/alvaradocx/Documents/single_cell/ALS_sumstats38.txt"
ALS_vg_gwas = "/Users/alvaradocx/Documents/single_cell/GCST90027164_buildGRCh37.tsv.gz"

# munge
AD_formatted <- MungeSumstats::format_sumstats(path=AD_gwas,
                                           save_path = "/Users/alvaradocx/Documents/single_cell/AD_Bellenguez.formatted.tsv",
                                           ref_genome ="GRCh38",
                                           force_new = TRUE)

PD_formatted <- format_sumstats(path=PD_gwas,  # for PD sums stats we want to liftover to GR38 since the other two dx sum stats +single cell expression data 
                                save_path = "/Users/alvaradocx/Documents/single_cell/PD_Nalls_no23.formatted.tsv", 
                                ref_genome ="GRCh37",
                                convert_ref_genome = "GRCh38",
                                force_new = TRUE)

PD_gwas_formatted <- format_sumstats(path=PD_meta_gwas, 
                                save_path = "/Users/alvaradocx/Documents/single_cell/PD_Nalls.formatted.tsv", 
                                ref_genome ="GRCh38",
                                force_new = TRUE)

LBD_formatted <- format_sumstats(path=LBD_gwas, 
                                 save_path = "/Users/alvaradocx/Documents/single_cell/LBD_Chia.formatted.tsv", 
                                 ref_genome ="GRCh38")

FTLD_formatted <- format_sumstats(path=FTLD_gwas, 
                                 save_path = "/Users/alvaradocx/Documents/single_cell/FTLD_Pottier.formatted.tsv", 
                                 ref_genome = "GRCh37",
                                 convert_ref_genome = "GRCh38")

PSP_formatted <- format_sumstats(path=PSP_gwas, 
                                 save_path = "/Users/alvaradocx/Documents/single_cell/PSP_Hoglinger.formatted.tsv", 
                                 ref_genome = "GRCh37",
                                 convert_ref_genome = "GRCh38",
                                 force_new = TRUE)
ALS_formatted <- format_sumstats(path=ALS_gwas, 
                                 save_path = "/Users/alvaradocx/Documents/single_cell/ALS_Nicolas.formatted.tsv", 
                                 ref_genome = "GRCh37",
                                 convert_ref_genome = "GRCh38")

ALS_vg_formatted <- format_sumstats(path=ALS_vg_gwas, 
                                 save_path = "/Users/alvaradocx/Documents/single_cell/ALS_vanRheenen.formatted.tsv", 
                                 ref_genome = "GRCh37",
                                 convert_ref_genome = "GRCh38")

ALS_vg
FTLD2_formatted <- format_sumstats(path=FTLD_gwas, 
                                  save_path = "/Users/alvaradocx/Documents/single_cell/FTLD2_Pottier.formatted.tsv", 
                                  ref_genome = "GRCh38")


# SNP to Gene Mapping
AD_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = AD_formatted,
  genome_build = "GRCh38",
  N = 788989)

LBD_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = LBD_formatted,
  genome_build = "GRCh38",
  N = 6618)

PD_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = PD_formatted,
  genome_build = "GRCh38",
  N = 1456306)

PD_gwas_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = PD_gwas_formatted,
  genome_build = "GRCh38",
  N = 1530403,
  force_new = TRUE)

FTLD_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = FTLD_formatted,
  genome_build = "GRCh38",
  N = 1345,
  force_new = TRUE)

ALS_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = ALS_formatted,
  genome_build = "GRCh38",
  N = 80610,
  force_new = TRUE)

ALS_vg_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = ALS_vg_formatted,
  genome_build = "GRCh38",
  N = 152268,
  force_new = TRUE)

PSP_gene <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = PSP_formatted,
  genome_build = "GRCh38",
  N = 8972,
  force_new = TRUE)


# alternative ALs GWAS from package
metagwas <- MungeSumstats::find_sumstats(traits = c("Amyotrophic lateral sclerosis"))
datasets <- MungeSumstats::import_sumstats(ids = "ebi-a-GCST005647",
                                           ref_genome = "GRCH37")
#save ALS GWAS from the ieu open GWAS project to a temp directory
ALSvcfPth <- "/var/folders/9f/hk0dc7qs6gq4v044z86h86vx5t281n/T//Rtmp0kaRkU/ebi-a-GCST005647/ebi-a-GCST005647.tsv.gz"
reformatted_vcf <- 
  MungeSumstats::format_sumstats(path=ALSvcfPth, 
                                 save_path = "/Users/alvaradocx/Documents/single_cell/ALS_Nicolas_mungess.formatted.tsv", 
                                 ref_genome="GRCh37",
                                 convert_ref_genome = "GRCh38",
                                 force_new = TRUE)

ALS_gene2 <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted = reformatted_vcf,
  genome_build = "GRCh38",
  N = 80610,
  force_new = TRUE)
