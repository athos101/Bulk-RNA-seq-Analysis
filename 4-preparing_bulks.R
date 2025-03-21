library("tidyverse")
library("purrr")
library("biomaRt")
library("org.Hs.eg.db")

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
readCount <- readRDS(file = "bulk_data/merged_bulks_clean.rds")