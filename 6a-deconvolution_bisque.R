## ABOUT THE FILE
#//=======================================================================//

# Author: Athos B. Schuck

## STEP 1: LOAD LIBRARIES
#//=======================================================================//
library(BisqueRNA)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(readr)
library(SingleCellExperiment)
library(DeconvoBuddies)

## STEP 2: LOAD DATA
#//=======================================================================//
# LUAD
bulkdata <- readRDS('bulk_data/merged_bulks_clean.rds')
bulkdata <- as.matrix(bulkdata)
bulk.eset <- ExpressionSet(assayData = bulkdata)
rm(bulkdata)
sobj <- readRDS(file = '/mnt/hd6tb/athos_hub/rds_data/atlas_clean.rds')

## STEP 3: FILTERING MARKERS WITH MEANRATIO
#//=======================================================================//
# First, we clean and prepare the data.
sobj <- SetIdent(sobj, value=sobj$celltype)
sobj <- subset(sobj, idents=c("Remove"), invert=TRUE)

sce <- as.SingleCellExperiment(sobj)
colData(sce) <- DataFrame(sobj@meta.data)
rowData(sce) <- rownames(sobj)
# ONLY IF THE MATRIX IS RAW!
assay(sce, "logcounts") <- log1p(assay(sce, "counts"))
# IF IT IS ALREADY NORMALIZED
assay(sce, "logcounts") <- assay(sce,"counts")
# THEN WE DO THE MEAN RATIO

marker_stats_MeanRatio <- get_mean_ratio(
  sce = sce, # sce is the SingleCellExperiment with our data
  assay_name = "counts", ## assay to use, we recommend logcounts [default]
  cellType_col = "celltype", # column in colData with cell type info
)

marker_stats_MeanRatio |>
  filter(MeanRatio.rank == 1)

# AND FINALLY, WE VERIFY THE GENES SELECTED

plot_gene_express(
  sce = sce,
  category = "celltype",
  genes = c("KLF6","CD52")
)

plot_marker_express(
  sce = sce,
  stats = marker_stats_MeanRatio,
  cell_type = "TCD4",
  n_genes = 10,
  cellType_col = "celltype"
)

## STEP 4: FILTERING BOTH BULK AND SINGLE-CELL WITH THE MARKERS
#//=======================================================================//
rownames(bulk.eset) <- rowData(bulk.eset)$Symbol
marker_genes <- marker_stats |>
  filter(MeanRatio.rank <= 25, gene %in% rownames(bulk.eset)) |>
  pull(gene)
sc_counts <- GetAssayData(sobj[marker_genes,], layer = "counts")
cellmeta <- sobj@meta.data
celltype <- cellmeta[,c("celltype", "patient")]
colnames(celltype) <- c("cellType", "SubjectName")
df <- data.frame(labelDescription=colnames(celltype), row.names = colnames(celltype))
#rm(sobj)
sc_annot <- new("AnnotatedDataFrame", data=celltype, varMetadata=df)
sc.eset <- ExpressionSet(assayData = as.matrix(sc_counts), phenoData = sc_annot)

#//=======================================================================//

res <- BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers=NULL, use.overlap=FALSE)
bisque_results <- as.data.frame(res$bulk.props)
bisque_results <- t(bisque_results)

#//=======================================================================//





































