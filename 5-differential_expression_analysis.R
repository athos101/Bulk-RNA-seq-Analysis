library("tidyverse")
library("dplyr")
library("biomaRt")

map_status <- function(status){
  if(status=="ADJ"){
    return("Adjacent")
  } else if (status=="TUM"){
    return("Tumor")
  }else if (status=="LINFO"){
    return("Lymph Node")
  }
}

process_count_data <- function(path_to_tabular, mart, status) {
  count_data <- read.csv(path_to_tabular, sep = "\t")
  file_name_info <- strsplit(basename(path_to_tabular),".")
  colnames(count_data)[2] <- paste("Pul_",file_name_info[1])
  file_name_info <- strsplit(basename(file_name_info[1]),"_")
  count_data$Status <- status
  count_data$Tissue <- map_status(file_name_info[2])
  count_data <- count_data[, sort(colnames(count_data))]
  return(count_data)
}


count_data$Geneid <- gsub("\\..*", "", count_data$Ensamble_id)
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = count_data$Ensamble_id,
  mart = mart
)
count_data <- merge(count_data, gene_annotations, by.x = "Ensamble_id", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(count_data)[colnames(count_data) == "external_gene_name"] <- "Gene"


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



