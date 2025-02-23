library("tidyverse")
library("dplyr")
library("purrr")
library("biomaRt")

# Lista de valores a serem buscados
valores <- duplicated_symbols

# Itera sobre cada valor e imprime os índices onde X tem esse valor
for (valor in valores) {
  indices <- rownames(count_data)[count_data$Geneid == valor]  # Obtém os índices
  print(paste0(valor, ": [", paste(indices, collapse = ", "), "]"))
}

ensambl_to_gene <- function(count_data, mart){
  count_data$Geneid <- gsub("\\..*", "", rownames(count_data))
  duplicated_ids <- any(duplicated(count_data$Geneid))
  
  if(duplicated_ids){
    cat("There are duplicated gene symbols in the data.\n")
    # Show the duplicated gene symbols
    duplicated_symbols <- count_data$Geneid[duplicated(count_data$Geneid)]
    print(duplicated_symbols)
    count_data <- count_data %>%
      group_by(Geneid) %>% 
      summarise(across(everything(), sum))
  } else {
    cat("There are no duplicated gene symbols in the data.\n")
  }
  
  count_data <- as.data.frame(count_data)
  rownames(count_data) <- count_data$Geneid
  
  # RETRIEVE BIOMART DATABASE INFORMATION
  gene_annotations <- getBM(
    # WHAT DATA DO WE WANT TO RETRIEVE
    attributes = c("ensembl_gene_id", "external_gene_name"),
    # FROM WHAT DATA WE WILL MAP FROM? ENSAMBL TO GENE SYMBOL!
    filters = "ensembl_gene_id",
    # WHERE ARE THE DATA THAT WE WILL MAP FROM
    values = count_data$Geneid,
    # DATABASE
    mart = mart
  )
  count_data <- merge(count_data, gene_annotations, by.x = "Geneid", by.y = "ensembl_gene_id", all.x = TRUE)
  count_data <- subset(count_data, select=-c(Geneid))
  colnames(count_data)[colnames(count_data) == "external_gene_name"] <- "Geneid"
  return(count_data)
}

## ACTUAL LOOP HERE

csv_data_folder <- "bulk_data"
target_path <- file.path(getwd(), csv_data_folder)

if(!dir.exists(target_path)){
  stop(paste("Folder not found: ", target_path))
}

condition_folders<-list.dirs(target_path, recursive=FALSE)

# MART - HUMAN GENES FOR ANNOTATION
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org") 

merged_counts <- data.frame()

for(subfolder in condition_folders){
  condition <- basename(subfolder)
  files_in_subfolder <- lapply(subfolder, list.files, full.names=TRUE)
  for(file in files_in_subfolder[[1]]){
    # READS THE FILE
    print(paste("Reading file: ",basename(file)))
    count_data <- read.csv(file, sep = "\t")
    # RENAME THE COUNTS COLUM TO THE NAME OF THE SAMPLE
    file_name_info <- strsplit(basename(file),"\\.")
    file_name_info <- paste("Pul_",file_name_info[[1]][1],"_",condition)
    colnames(count_data)[2] <- file_name_info
    rownames(count_data)<-count_data$Geneid
    count_data<-subset(count_data, select = -c(Geneid))
    if(dim(merged_counts)[1]==0){
      merged_counts<-count_data
    }else{
      merged_counts<-cbind(merged_counts, count_data)
    }
  }
}

merged_counts <- ensambl_to_gene(merged_counts, mart)

## CODE OF WILCOXON + EDGER TMM
library(edgeR)

# read data
readCount <- read.table(file = "examples/examples.countMatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
conditions <- read.table(file = "examples/examples.conditions.tsv", header = F)
conditions <- factor(t(conditions))
# edgeR TMM normalize
y <- DGEList(counts = readCount, group = conditions)
## Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep, keep.lib.sizes = FALSE]
## Perform TMM normalization and convert to CPM (Counts Per Million)
y <- calcNormFactors(y, method = "TMM")
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)
# Run the Wilcoxon rank-sum test for each gene
pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)
})
fdr <- p.adjust(pvalues, method = "fdr")
# Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
write.table(outRst[outRst$FDR<fdrThres,], file = "examples/examples.WilcoxonTest.rst.tsv", sep="\t", quote = F, row.names = T, col.names = T)


