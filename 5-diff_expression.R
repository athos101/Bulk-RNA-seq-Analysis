library("tidyverse")
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

build_conditions <- function(names){
  conditions <- data.frame(1)
  for(name in names){
    conditions[name] <- strsplit(name, "_")[[1]][3]
  }
  conditions <- conditions[2:(length(names)+1)]
  conditions <- factor(t(conditions))
  return(conditions)
}

# //-=-=-=-=-=-=-=-=-PREPARING ADJ SAMPLES FOR WILCX =-=-=-=-=-=-=-=
adj_samples <- lapply(colnames(readCount), FUN = function(name){
  if(length(grep("ADJ", name)))
    return(name)
})
adj_samples <- Filter(Negate(is.null), adj_samples)
bulk_adjs <- readCount[, (colnames(readCount) %in% adj_samples)]

conditions_adj <- build_conditions(colnames(bulk_adjs))

# //-=-=-=-=-=-=-=-=-PREPARING TUM SAMPLES FOR WILCX =-=-=-=-=-=-=-=

tum_samples <- lapply(colnames(readCount), FUN = function(name){
  if(length(grep("TUM", name)))
    return(name)
})

tum_samples <- Filter(Negate(is.null), tum_samples)
bulk_tums <- readCount[, (colnames(readCount) %in% tum_samples)]

conditions_tum <- build_conditions(colnames(bulk_tums))

# //-=/-=/-=/-=/-=/-=/ WILCOXON RANK SUM \=-\=-\=-\=-\=-\=-\=-\=-||

# Here, we create the conditions for our Wilcoxom runk sum,
# which contains for each sample the conditions to be compared.

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
    p <- wilcox.test(gene~conditions, data, exact=FALSE, correct=TRUE)$p.value
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
  print(min(outRst$FDR, na.rm=TRUE))
  fdrThres <- 0.05
  print(head(outRst[outRst$FDR<fdrThres,]))
  return(outRst[outRst$FDR<fdrThres,])
}

wilcox_adj <- generate_wilcoxon_table(bulk_adjs, conditions_adj)
write.table(wilcox_adj, file = "WCX_adj_rst.tsv", sep="\t", quote = F, row.names = T, col.names = T)
wilcox_tum <- generate_wilcoxon_table(bulk_tums, conditions_tum)
write.table(wilcox_tum, file = "WCX_tum_rst.tsv", sep="\t", quote = F, row.names = T, col.names = T)


# //-=/-=/-=/-=/-=/-=/ DESEQ2 \=-\=-\=-\=-\=-\=-\=-\=-||
library("DESeq2")
# identify biological replicates
deseq_conditions <- lapply(colnames(readCount), function(name){
  pieces <- str_split(name, "_")
  name <- paste(pieces[[1]][2], pieces[[1]][3], sep="_")
  return(name)
})

# assign replicates to each sample name to construct colData
my_colData <- as.data.frame(deseq_conditions)
my_colData <- t(my_colData)
rownames(my_colData) <- colnames(readCount)
colnames(my_colData)[1] <- "condition"

dds <- DESeqDataSetFromMatrix(countData = readCount,
                              colData = my_colData,
                              design = ~deseq_conditions)

dds <- DESeq(dds)

# //-=/-=/-=/-=/-=/-=/  GSVA  \=-\=-\=-\=-\=-\=-\=-\=-||
library(GSVA)

# LOAD FUNCTIONAL GENE EXPRESSION SIGNATURES
fges <- readRDS('E:/aula tumores/FGES_list.rds')
fges <- fges[1:27]

#Perform log-CPM normalization with edgeR
bulkdata <- cpm(as.matrix(counts_matrix_cleaned2), log = T)
#Convert matrix to ExpressionSet
bulkdata <- ExpressionSet(assayData = bulkdata)

#Run GSVA
#Gaussian for log CPM normalization, if raw counts provided used "Poisson"
bulkPar <- gsvaParam(bulkdata, fges, maxDiff=FALSE, kcdf = "Gaussian")
bulk_es <- gsva(bulkPar)

#Plot Heatmap
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]
pheatmap(exprs(bulk_es), col = hmcol, scale = "row",
         cluster_cols = T, cluster_rows = T,
         clustering_method = "ward.D2", angle_col = "45")