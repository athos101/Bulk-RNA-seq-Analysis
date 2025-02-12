library("DESeq2")
library("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

data <- read.csv("raw_counts.csv", header = T, row.names = "ensembl_id")
data <- data[,sort(colnames(data))]
head(data)

# TUTORIAL - VIDEO 

count_data <- read.csv("bulk_data/auto_anti_pdl1/58_ADJ.tabular", sep = "\t")
colnames(count_data)[2] <- "counts"
count_data$condition <- "Anti-PDL1"
count_data$tissue <- "Adjacent"
head(count_data)

count_data2 <- read.csv("bulk_data/auto_anti_pdl1/58_TUM.tabular", sep = "\t")
colnames(count_data2)[2] <- "counts"
count_data2$condition <- "Anti-PDL1"
count_data2$tissue <- "Tumor"
head(count_data2)

library(EnsDb.Hsapiens.v79)

# 1. Convert from ensembl.gene to gene.symbol
ensembl.genes <- c("ENSG00000150676", "ENSG00000099308", "ENSG00000142676", "ENSG00000180776", "ENSG00000108848", "ENSG00000277370", "ENSG00000103811", "ENSG00000101473")

geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

result_table

# Create dataframe
df <- data.frame(id=c(1,2,3,NA),
                 address=c('Orange St','Anton Blvd','Jefferson Pkwy',''),
                 work_address=c('Main St',NA,'Apple Blvd','Portola Pkwy'))
df$address[df$address == 'Orange St'] <- 'Portola Pkwy'



