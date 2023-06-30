# library loads
library(EWCE)
library(sctransform)
library(readxl)
library(stats)
library(tidyverse)


# load in files
init_matrix <- read.csv("/Users/alvaradocx/Documents/single_cell/genexcluster_matrix_nodups.csv", sep = ',', row.names = 1) # read in expression matrix
colnames(init_matrix)[1]<- "HGNC" # change gene column name
gene_list <- init_matrix$HGNC
rownames(init_matrix) <- init_matrix$HGNC  # Make hgnc column row names
init_matrix <- init_matrix[c(-1)]

cluster_meta <- read.csv("/Users/alvaradocx/Documents/single_cell/cluster_metadata.csv") # cluster metadata
cluster_meta$Cluster <- sub("^", "V", cluster_meta$Cluster) # add "V" in front of each cluster for ease of processing
cluster_meta <- cluster_meta[-c(462),] # remove last row since it is empty
# Replace na values with string placeholder otherwise drop_uninformative_genes function will throw an error
cluster_meta$`Class auto-annotation` <- cluster_meta$`Class auto-annotation` %>% replace_na('none')

# create main list object
main <- list()

# create matrix object from expression matrix
mod_matrix <- as.matrix(init_matrix)

# Correct gene symbols
main$exp = fix_bad_hgnc_symbols(mod_matrix) # add corrected gene symbols output to main list object

# add metadata to main list object
main$annot <- cluster_meta

# Calculate specificity matrices
exp_DROPPED <- drop_uninformative_genes(
  exp = main$exp, 
  input_species = "human",
  output_species = "human",
  level2annot = main$annot$Class.auto.annotation) 

annotLevels <- list(level1class = main$annot$Supercluster,
                    level2class = main$annot$Class.auto.annotation)

# Generate CellTypeDataset object
adult_ctd <- generate_celltype_data(
  exp = exp_DROPPED,
  annotLevels = annotLevels,
  groupName = "adult_cluster_genes_mod",
  savePath = "/Users/alvaradocx/Documents/single_cell",
  return_ctd=TRUE) 


# bootstrap enrichment
boot <- bootstrap_enrichment_test(sct_data = adult_ctd, hits = gene_list, reps = 10000, annotLevel = 1, sctSpecies = 'human', output_species = 'human')
