library("tidyverse")
library("dplyr")
library("purrr")
library("biomaRt")

ensambl_to_gene <- function(count_data, mart){
  count_data$Geneid <- gsub("\\..*", "", count_data$Geneid)
  duplicated_ids <- any(duplicated(count_data$Geneid))

  if(duplicated_ids){
    count_data <- count_data %>%
      group_by(Geneid) %>% 
      summarise(across(everything(), sum))
  }
  
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
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

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
    colnames(count_data)[2] <- paste("Pul_",file_name_info[[1]][1])
    rownames(count_data)<-count_data$Geneid
    # WE SORT THE COLUMNS TO FURTHER MERGE THEM
    count_data <- count_data[, sort(colnames(count_data))]
    if(dim(merged_counts)[1]==0){
      merged_counts<-count_data
    }else{
      merged_counts<-cbind(merged_counts, count_data)
    }
  }
}

merged_counts <- ensambl_to_gene(count_data, mart)



