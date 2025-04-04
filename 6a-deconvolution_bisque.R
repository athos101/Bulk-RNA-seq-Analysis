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
library(Matrix)
library(SingleCellExperiment)
library(DeconvoBuddies)

## STEP 2: LOAD DATA
#//=======================================================================//
# LUAD
bulkdata <- readRDS('bulk_data/merged_bulks_clean.rds')
bulkdata <- as.matrix(bulkdata)
bulk.eset <- ExpressionSet(assayData = bulkdata)
rm(bulkdata)
sobj <- readRDS(file = '/mnt/hd6tb/athos_hub/rds_data/atlas_raw_18k.rds')

## STEP 3: FILTERING MARKERS WITH MEANRATIO
#//=======================================================================//
# First, we clean and prepare the data.
sobj <- SetIdent(sobj, value=sobj$celltype)
sobj <- subset(sobj, idents=c("Remove"), invert=TRUE)
sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

# We extract only the first 5000 most variable genes for computational efficiency.
sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 5000)
sobj <- subset(sobj, features = VariableFeatures(sobj))
gc()

sce <- as.SingleCellExperiment(sobj)
colData(sce) <- DataFrame(sobj@meta.data)
rowData(sce) <- rownames(sobj)
assay(sce, "logcounts") <- log1p(assay(sce, "counts"))
# THEN WE DO THE MEAN RATIO

marker_stats_MeanRatio <- get_mean_ratio(
  sce = sce, # sce is the SingleCellExperiment with our data
  assay_name = "logcounts", ## assay to use, we recommend logcounts [default]
  cellType_col = "celltype", # column in colData with cell type info
)

marker_stats_MeanRatio

marker_stats_MeanRatio |>
  filter(MeanRatio.rank == 1)

# AND FINALLY, WE VERIFY THE GENES SELECTED

plot_gene_express(
  sce = sce,
  category = "cellType_broad_hc",
  genes = c("KRT19", "MK167")
)

plot_marker_express(
  sce = sce,
  stats = marker_stats_MeanRatio,
  cell_type = "Excit",
  n_genes = 10,
  cellType_col = "cellType.target"
)

## STEP 4: CLEANING THE SC RNA SEQ DATA
#//=======================================================================//

sobj <- sobj[, colData(sobj)$cell_type != "Remove"]
ct_tab <- table(sobj$cell_type)
my_cell_types <- names(ct_tab[ct_tab > 50])
sobj <- sobj[,sobj$cell_type %in% my_cell_type]

## STEP 5: OBTAINING THE REPRESENTATIVE GENES OF EACH CELL TYPE
#//=======================================================================//
gene_sums <- rowSums(logcounts(sobj))
top_gene <- gene_sums > median(gene_sums)
sobj <- sobj[top_gene,]
marker_stats <- get_mean_ratio2(sobj, cellType_col = "cellType_broad")
marker_stats

DeconvoBuddies::plot_marker_express(sobj,
                                    stats = marker_stats, 
                                    cell_type = 'MYCELLTYPE', 
                                    cellType_col = "cell_type", 
                                    n_genes = 10, 
                                    rank_col = "rank_ratio",
                                    anno_col = "anno_ratio",
)

DeconvoBuddies::plot_gene_express(sobj,
                                  genes = c("MYGENE1","MYGENE2"))

## STEP 6: PREPARING SINGLE-CELL COUNTS
#//=======================================================================//

sc_counts <- GetAssayData(sobj, layer = "counts")
cellmeta <- sobj@meta.data
celltype <- cellmeta[,c("celltype", "patient")]
colnames(celltype) <- c("cellType", "SubjectName")
df <- data.frame(labelDescription=colnames(celltype), row.names = colnames(celltype))
#rm(sobj)
sc_annot <- new("AnnotatedDataFrame", data=celltype, varMetadata=df)
sc.eset <- ExpressionSet(assayData = as.matrix(sc_counts), phenoData = sc_annot)

## STEP 7: PREPARING SINGLE-CELL COUNTS
#//=======================================================================//



git commit -m "minor changes"


# A tibble: 16 × 8
gene     cellType.target mean.target cellType.2nd  mean.2nd MeanRatio MeanRatio.rank MeanRatio.anno                    
<chr>    <chr>                 <dbl> <chr>            <dbl>     <dbl>          <int> <chr>                             
  1 KRT19    Epithelial            1.64  Fibroblast      0.195      8.43               1 Epithelial/Fibroblast: 8.426      
2 MKI67    T cycle               1.39  Epithelial      0.178      7.80               1 T cycle/Epithelial: 7.795         
3 CD69     TCD4                  1.17  T naive         1.19       0.987              1 TCD4/T naive: 0.987               
4 IL2RA    Treg                  0.726 T cycle         0.338      2.15               1 Treg/T cycle: 2.147               
5 CCL5     TCD8                  2.28  NK              1.97       1.16               1 TCD8/NK: 1.155                    
6 LEPROTL1 T naive               0.870 T cycle         0.723      1.20               1 T naive/T cycle: 1.204            
7 KLRF1    NK                    0.752 TCD8            0.0645    11.7                1 NK/TCD8: 11.67                    
8 SFRP2    Fibroblast            1.78  Smooth Muscle   0.0879    20.3                1 Fibroblast/Smooth Muscle: 20.261  
9 MS4A1    B cells               1.64  Endothelial     0.0579    28.4                1 B cells/Endothelial: 28.405       
10 MS4A4A   Myeloid cells         1.38  DCs             0.269      5.14               1 Myeloid cells/DCs: 5.143          
11 CD1C     DCs                   1.17  B cells         0.115     10.1                1 DCs/B cells: 10.149               
12 G0S2     Neutrophils           2.68  DCs             1.10       2.44               1 Neutrophils/DCs: 2.435            
13 CLEC4C   pDC                   0.819 DCs             0.0119    68.7                1 pDC/DCs: 68.671                   
14 EMCN     Endothelial           1.30  Smooth Muscle   0.0113   115.                 1 Endothelial/Smooth Muscle: 115.133
15 IGHGP    Plasma cells          1.97  B cells         0.0479    41.2                1 Plasma cells/B cells: 41.163      
16 RGS5     Smooth Muscle         2.21  Endothelial     0.374      5.90               1 Smooth Muscle/Endothelial: 5.898  




# A tibble: 1,089 × 8
gene   cellType.target mean.target cellType.2nd  mean.2nd MeanRatio MeanRatio.rank MeanRatio.anno                 
<chr>  <chr>                 <dbl> <chr>            <dbl>     <dbl>          <int> <chr>                          
  1 KRT19  Epithelial            1.64  Fibroblast       0.195     8.43               1 Epithelial/Fibroblast: 8.426   
2 KRT8   Epithelial            1.38  Fibroblast       0.211     6.57               2 Epithelial/Fibroblast: 6.571   
3 SLPI   Epithelial            1.37  Fibroblast       0.296     4.65               3 Epithelial/Fibroblast: 4.646   
4 CLDN7  Epithelial            0.875 Myeloid cells    0.231     3.79               4 Epithelial/Myeloid cells: 3.793
5 KRT18  Epithelial            1.33  Smooth Muscle    0.445     3.00               5 Epithelial/Smooth Muscle: 2.996
6 NQO1   Epithelial            0.956 Endothelial      0.319     3.00               6 Epithelial/Endothelial: 2.996  
7 EIF2S2 Epithelial            1.20  Fibroblast       0.995     1.20               7 Epithelial/Fibroblast: 1.204   
8 EFNA1  Epithelial            0.822 Endothelial      0.688     1.19               8 Epithelial/Endothelial: 1.193  
9 JPT1   Epithelial            0.943 T cycle          0.936     1.01               9 Epithelial/T cycle: 1.008      
10 HSPD1  Epithelial            1.33  Fibroblast       1.34      0.994             10 Epithelial/Fibroblast: 0.994   
















































