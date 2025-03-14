library("tidyverse")
library("dplyr")
library("edgeR")

readCount <- readRDS(file = "bulk_data/merged_bulks_clean.rds")

## REPRODUCING THE WILCOXON + EDGER TMM

get_condition <- function(string){
  if(grepl("ADJ",string, fixed = TRUE))
    return("Adjacent")
  if(grepl("TUM",string, fixed = TRUE))
    return("Tumor")
  if(grepl("LINFO",string, fixed = TRUE))
    return("LymphNode")
}

# Here, we create the conditions for our Wilcoxom runk sum,
# which contains for each sample the conditions to be compared.

adj_samples <- lapply(colnames(readCount), FUN = function(name){
  if(length(grep("ADJ", name)))
    return(name)
})
adj_samples <- Filter(Negate(is.null), adj_samples)
bulk_adjs <- readCount[, (colnames(readCount) %in% adj_samples)]

tum_samples <- lapply(colnames(readCount), FUN = function(name){
  if(length(grep("TUM", name)))
    return(name)
})
tum_samples <- Filter(Negate(is.null), tum_samples)
bulk_tums <- readCount[, (colnames(readCount) %in% tum_samples)]

linfo_samples <- lapply(colnames(readCount), FUN = function(name){
  if(length(grep("LINFO", name)))
    return(name)
})
linfo_samples <- Filter(Negate(is.null), linfo_samples)
bulk_linfo <- readCount[, (colnames(readCount) %in% linfo_samples)]

conditions <- c("")
conditions <- factor(t(conditions))

## CODE OF WILCOXON + EDGER TMM
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
  return(outRst[outRst$FDR<fdrThres])
}