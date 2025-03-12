library("tidyverse")
library("dplyr")
library("purrr")
library("biomaRt")
library("org.Hs.eg.db")
library("edgeR")

remove_ensambl_suffix <- function(count_data){
  count_data$Geneid <- gsub("\\..*", "", rownames(count_data))
  duplicated_ids <- any(duplicated(count_data$Geneid))
  
  if(duplicated_ids){
    count_data <- count_data %>%
      group_by(Geneid) %>% 
      summarise(across(everything(), sum))
      remove(duplicated_ids)
  }
  
  count_data <- as.data.frame(count_data)
  rownames(count_data) <- count_data$Geneid
  count_data <- subset(count_data, select=-c(Geneid))
  return(count_data)
}

ensambl_to_gene <- function(count_data){
  count_data <- remove_ensambl_suffix(count_data)
  count_data$symbol <- mapIds(org.Hs.eg.db, keys = rownames(count_data), keytype = "ENSEMBL", column = "SYMBOL")
  filtered_counts <- count_data
  na_count <- sum(is.na(filtered_counts))
  cat("Number of NAs in the gene symbols column:", na_count, "\n")
  filtered_counts <- filtered_counts[!is.na(filtered_counts$symbol),]
  dim(filtered_counts)
  
  duplicated_genes <- any(duplicated(filtered_counts$symbol))
  if(duplicated_genes){
    cat("There are duplicated gene symbols in the data. \n")
    duplicated_symbols <- filtered_counts$symbol[duplicated(filtered_counts$symbol)]
    print(duplicated_symbols)
  } else {
    cat("There are no duplicated gene symbols in the data.\n")
  }
  
  duplicates <- duplicated(filtered_counts$symbol) | duplicated(filtered_counts$symbol, fromLast = TRUE)
  
  filtered_counts <- filtered_counts %>%
    group_by(symbol) %>%
    summarise(across(everything(), sum))
  
  count_data <- as.data.frame(filtered_counts)
  return(count_data)
}

ensambl_to_gene_alt <- function(count_data, mart){
  count_data<-remove_ensambl_duplicates(count_data)
  
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
csv_data_folder <- file.path(getwd(), csv_data_folder)

if(!dir.exists(csv_data_folder)){
  stop(paste("Folder not found: ", csv_data_folder))
}

condition_folders <- list.dirs(csv_data_folder, recursive=FALSE)

bulk_list <- list()

for(subfolder in condition_folders){
  condition <- basename(subfolder)
  csv_files <- list.files(path = subfolder, pattern = "\\.tabular$", full.names = TRUE)
  
  for(file in csv_files){
    # READS THE FILE
    print(paste("Reading file: ",basename(file)))
    count_data <- read.csv(file, sep = "\t")
    # RENAME THE COUNTS COLUM TO THE NAME OF THE SAMPLE
    file_name_info <- strsplit(basename(file),"\\.")
    file_name_info <- paste("Pul",file_name_info[[1]][1],"_",condition, sep = "")
    colnames(count_data)[2] <- file_name_info
    rownames(count_data)<-count_data$Geneid
    count_data<-subset(count_data, select = -c(Geneid))
    
    if(length(grep("ENSG", rownames(count_data)[1],fixed = TRUE))){
      # MART - HUMAN GENES FOR ANNOTATION - only for the ensambl_to_gene_alt
      #mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
      count_data <- ensambl_to_gene(count_data = count_data)
    } else {
      count_data$symbol <- rownames(count_data)
      count_data <- as.data.frame(count_data)
    }
    
    bulk_list[[file_name_info]]<-count_data
  }
}

merged_bulks <- Reduce(function(x,y) merge(x,y,by="symbol", all=TRUE), bulk_list)
rownames(merged_bulks) <- merged_bulks$symbol
merged_bulks <- subset(merged_bulks, select=-c(symbol))

# Some ensambl_ids have duplicated values after removing the suffix.
# It was found that all these suffixes were PARY genes, homologous genes
# of chromossomes X and Y. We just sum them. You don't need to execute this
# on the pipeline (:
valores <- duplicated_symbols
for (valor in valores) {
  indices <- rownames(count_data)[count_data$Geneid == valor]  # Obtém os índices
  print(paste0(valor, ": [", paste(indices, collapse = ", "), "]"))
}

merged_bulks[is.na(merged_bulks)] <- 0

saveRDS(merged_bulks, file = 'bulk_data/merged_bulks_clean.rds')
readCount <- readRDS(file = "merged_bulks_clean.rds")

## REPRODUCING THE WILCOXON + EDGER TMM

get_condition <- function(string){
  if(grepl("ADJ",string, fixed = TRUE))
    return("Adjacent")
  if(grepl("TUM",string, fixed = TRUE))
    return("Tumor")
  if(grepl("LINFO",string, fixed = TRUE))
    return("LymphNode")
}

conditions <- data.frame(1)
for(col in colnames(readCount)){
  #conditions[col] <- get_condition(col)
  conditions[col] <- str_split(col, "_")[[1]][3]
}
conditions <- conditions[2:(length(colnames(readCount))+1)]
conditions <- factor(t(conditions))

## CODE OF WILCOXON + EDGER TMM

# Read count matrix: genes for rows and samples for clumns.
readCount <- read.table(file = "examples/examples.countMatrix.tsv", header = T, row.names = 1, stringsAsFactors = F, check.names = F)
# Conditions is a single row with the condition labels for each sample present in the read count matrix.
conditions <- read.table(file = "examples/examples.conditions.tsv", header = F)
conditions <- factor(t(conditions))

generate_wilcoxon_table <- function(readCount, conditions){
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
  
  # IF >0, IT IS MORE EXPRESSED IN DATACON2. ELSE, IT IS MORE EXPRESSED IN DATACON1. 
  foldChanges <- log2(rowMeans(dataCon1)/rowMeans(dataCon2))
  # Output results based on the FDR threshold 0.05
  outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
  rownames(outRst) <- rownames(count_norm)
  outRst <- na.omit(outRst)
  fdrThres <- 0.05
  write.table(outRst[outRst$FDR<fdrThres,], file = "WilcoxonTest_Clinical.rst.tsv", sep="\t", quote = F, row.names = T, col.names = T)
}


