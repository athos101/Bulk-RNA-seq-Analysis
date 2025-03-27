## ABOUT THE FILE
#//=======================================================================//

# Author: Athos B. Schuck

## STEP 1: LOAD LIBRARIES
#//=======================================================================//
library(BisqueRNA)
library(SeuratObject)
library(ggplot2)
library(DeconvoBuddies)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(readr)
library(Matrix)
library(DeconvoBuddies)

## STEP 2: LOAD DATA
#//=======================================================================//

bulkdata <- readRDS('bulk_data/merged_bulks_clean.rds')
bulkdata <- as.matrix(bulkdata)
bulk.eset <- ExpressionSet(assayData = bulkdata)
rm(bulkdata)
sobj <- readRDS(file = '../R_Datasets/LUSC_subset_raw.rds')


## STEP 3: FILTERING MARKERS WITH MEANRATIO
#//=======================================================================//

# FIRST WE FILTER THE NUMBER OF GENES TO COMPUTATIONAL EFFICIENCY
markers_df <- read.csv('/home/usuario/Athos Plaza/markers_LUAD_L3_top50.csv')
rownames(markers_df) <- markers_df$X
markers_df <- markers_df[ , !names(markers_df) == "X"]

gene_list <- unlist(markers_df)
gene_list <- unique(gene_list)
length(gene_list)

genes_present <- gene_list[gene_list %in% rownames(sobj)]
genes_not_present <- setdiff(gene_list, genes_present)
length(genes_present)

sobj <- subset(sobj, features = genes_present)
sobj <- as.SingleCellExperiment(sobj)

# THEN WE DO THE MEAN RATIO

marker_stats_MeanRatio <- get_mean_ratio(
  sce = sce, # sce is the SingleCellExperiment with our data
  assay_name = "logcounts", ## assay to use, we recommend logcounts [default]
  cellType_col = "annotation_X", # column in colData with cell type info
  gene_ensembl = NULL, # column in rowData with ensembl gene ids
  gene_name = "gene_name" # column in rowData with gene names/symbols
)

marker_stats_MeanRatio

marker_stats_MeanRatio |>
  filter(MeanRatio.rank == 1)

# AND FINALLY, WE VERIFY THE GENES SELECTED

plot_gene_express(
  sce = sce,
  category = "cellType_broad_hc",
  genes = c("MYGENE1", "MYGENE2")
)

plot_marker_express(
  sce = sce,
  stats = marker_stats,
  cell_type = "Excit",
  n_genes = 10,
  cellType_col = "cellType_broad_hc"
)


































































