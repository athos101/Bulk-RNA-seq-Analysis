library("DESeq2")
library("tidyverse")
library("org.Hs.eg.db")
library("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# READING THE CSV FILES
count_data <- read.csv("bulk_data/auto_anti_pdl1/58_ADJ.tabular", sep = "\t", row.names="Geneid")
# RENAMING COUNT COLUMN TO COUNTS
colnames(count_data)[1] <- "counts"
# ADDING CONDITION AND TISSUE INFO
count_data$condition <- "Anti-PDL1"
count_data$tissue <- "Adjacent"
# SORTING COLUMNS SO THAT ALL CSV COLUMNS READ WILL BE IN THE SAME ORDER
count_data <- count_data[,sort(colnames(count_data))]
head(count_data)

count_data2 <- read.csv("bulk_data/auto_anti_pdl1/58_TUM.tabular", sep = "\t", row.names="Geneid")
# RENAMING COUNT COLUMN TO COUNTS
colnames(count_data2)[1] <- "counts"
# ADDING CONDITION AND TISSUE INFO
count_data2$condition <- "Anti-PDL1"
count_data2$tissue <- "Adjacent"
# SORTING COLUMNS SO THAT ALL CSV COLUMNS READ WILL BE IN THE SAME ORDER
count_data2 <- count_data2[,sort(colnames(count_data2))]
head(count_data2)
# CONVERTING ENSAMBL TO GENE SYMBOL
if(substr(rownames(count_data2)[1],1,4)=="ENSG"){
  # CLEAN SUFFIX
  rownames(count_data2) <- gsub("\\..*", "", rownames(count_data2))
  head(count_data2)
  # CONVERT ENSEMBL TO SYMBOL
  count_data2 <- mapIds(org.Hs.eg.db,keys=count_data2$Ensembl,keytsub(ype = "ENSEMBL",column="SYMBOL"))
}

# 1. Convert from ensembl.gene to gene.symbol
colnames(count_data2)[1] <- "ensembl"
count_data2 <- mapIds(org.Hs.eg.db,keys=count_data2$Ensembl,keytsub(ype = "ENSEMBL",column="SYMBOL"))
?mapIds








