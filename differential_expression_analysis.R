library("DESeq2")
library("tidyverse")

data <- read.csv("raw_counts.csv", header = T, row.names = "ensembl_id")
data <- data[,sort(colnames(data))]
head(data)
colSums(data)

# identify biological replicates
condition <- c(rep("LNCAP_Hypoxia", 2), 
               rep("LNCAP_Normoxia", 2), 
               rep("PC3_Hypoxia", 2), 
               rep("PC3_Normoxia", 2))

my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(data)
my_colData

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = my_colData,
                              design = ~condition)

dds <- DESeq(dds)
dds

head(dds@assays@data$counts)

normalized_counts <- counts(dds, normalized = T)

head(normalized_counts)
annotation <- read.csv("GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
head(annotation)



