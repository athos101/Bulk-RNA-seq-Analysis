library("tidyverse")
library("dplyr")
library("biomaRt")

process_count_data <- function(path_to_tabular, mart) {
  count_data <- read.csv(path_to_tabular, sep = "\t")
  colnames(count_data)[2] <- "Counts"
  count_data$Status <- "Anti-PDL1"
  count_data$Tissue <- "Adjacent"
  count_data <- count_data[, sort(colnames(count_data))]
  if (substr(count_data$Geneid[1], 1, 4) == "ENSG") {
    count_data$Geneid <- gsub("\\..*", "", count_data$Geneid)
    gene_annotations <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
      filters = "ensembl_gene_id",
      values = count_data$Geneid,
      mart = mart
    )
    count_data <- merge(count_data, gene_annotations, by.x = "Geneid", by.y = "ensembl_gene_id", all.x = TRUE)
    colnames(count_data)[colnames(count_data) == "external_gene_name"] <- "Gene"
  }
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
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

for(subfolder in condition_folders){
  condition <- basename(subfolder)
  files_in_subfolder <- lapply(subfolder, list.files, full.names=TRUE)
  for(file in files_in_subfolder){
    count_data <- process_count_data(file, mart)
    sample_name <- basename(file)
  }
}



