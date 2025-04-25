library("tidyverse")
library("edgeR")
library("GSVA")
library("SingleCellExperiment")
library("RColorBrewer")

bulkdata <- readRDS(file = "bulk_data/merged_bulks_clean.rds")

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


# //-=/-=/-=/-=/-=/-=/-=/ DESEQ2 \=-\=-\=-\=-\=-\=-\=-\=-\=-\=-||
library("DESeq2")
# identify biological replicates
deseq_conditions <- lapply(colnames(bulkdata), function(name){
  pieces <- str_split(name, "_")
  name <- paste(pieces[[1]][2], pieces[[1]][3], sep="_")
  return(name)
})

# assign replicates to each sample name to construct colData
my_colData <- as.data.frame(deseq_conditions)
my_colData <- t(my_colData)
rownames(my_colData) <- colnames(bulkdata)
colnames(my_colData)[1] <- "deseq_conditions"
deseq_conditions <- unlist(deseq_conditions)
deseq_conditions <- as.factor(deseq_conditions)

dds <- DESeqDataSetFromMatrix(countData = bulkdata,
                              colData = my_colData,
                              design = ~deseq_conditions)

dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = T)
vsd <- vst(dds, blind = TRUE)
plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( vsd.obj$condition )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
}
plotDists(vsd)

variable_gene_heatmap <- function (vsd.obj, num_genes = 500, title = "") {
  brewer_palette <- "RdBu"
  # Ramp the color in order to get the scale.
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  # get the stabilized counts from the vsd object
  stabilized_counts <- assay(vsd.obj)
  # calculate the variances by row(gene) to find out which genes are the most variable across the samples.
  row_variances <- rowVars(stabilized_counts)
  # get the top most variable genes
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  # subtract out the means from each row, leaving the variances for each gene
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  # reconstruct colData without sizeFactors for heatmap labeling
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  # draw heatmap using pheatmap
  pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 8, fontsize_row = 300/num_genes, border_color = NA, main = title)
}

variable_gene_heatmap(vsd, num_genes = 40)

# //-=/-=/-=/-=/-=/-=/-=/ PLOTTING PCA \=-\=-\=-\=-\=-\=-\=-\=-\=-\=-||

library(ggrepel)
plot_PCA = function (vsd.obj) {
  pcaData <- plotPCA(vsd.obj,  intgroup = c("deseq_conditions"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=deseq_conditions)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}
plot_PCA(vsd)

# //-=/-=/-=/-=/-=/-=/ PLOTTING DEG IN THE SAME TISSUE \=-\=-\=-\=-\=-\=-\=-||

generate_DESeq_object <- function (my_data, groups) {
  data_subset1 <- my_data[,grep(groups[1], colnames(my_data))]
  data_subset2 <- my_data[,grep(groups[2], colnames(my_data))]
  my_countData <- cbind(data_subset1, data_subset2)
  condition <- c(rep(groups[1],ncol(data_subset1)), rep(groups[2],ncol(data_subset2)))
  my_colData <- as.data.frame(condition)
  rownames(my_colData) <- colnames(my_countData)
  print(my_colData)
  dds <- DESeqDataSetFromMatrix(countData = my_countData,
                                colData = my_colData,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = T)
  return(dds)
}

tum_samples <- generate_DESeq_object(readCount, c("TUM_Control", "TUM_aPDL1"))
adj_samples <- generate_DESeq_object(readCount, c("ADJ_Control", "ADJ_aPDL1"))

tum_samples <- vst(tum_samples, blind = T)
adj_samples <- vst(adj_samples, blind = T)

a <- variable_gene_heatmap(tum_samples, 30, title = "Tumor variable genes")

b <- variable_gene_heatmap(adj_samples, 30, title = "Adjacent variable genes")

gridExtra::grid.arrange(a[[4]],b[[4]], nrow = 1)

# //-=/-=/-=/-=/-=/-=/ GENE EXPRESSION BETWEEN SAMPLES \=-\=-\=-\=-\=-\=-\=-||
# Find rownames that start with "IG"
ig_genes <- rownames(dds)[grep("^IGHA", rownames(dds))]

# Print the gene names
print(ig_genes)

gene <- "IGHG4"

normalized_data <- counts(dds, normalized = T)
expression <- normalized_data[gene,]
#condition <- rownames(dds@colData)
condition <- dds@colData$deseq_conditions
gene_tib <- tibble(condition = condition, expression = expression)

ggplot(gene_tib, aes(x = condition, y = expression))+
  geom_boxplot(outlier.size = 0)+
  geom_point()+
  labs (title = paste0("Expression of ", gene), x = "group", y = paste0("Normalized expression (", "Deseq2" , ")"))+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12))

bulk_data <- readRDS(file = "bulk_data/merged_bulks_clean.rds")

# //-=/-=/-=/-=/-=/-=/  GSVA  \=-\=-\=-\=-\=-\=-\=-\=-||
library(GSVAdata)
library("pheatmap")
data(c2BroadSets)

# LOAD FUNCTIONAL GENE EXPRESSION SIGNATURES
fges <- read.csv('bulk_data/FGES_29.csv')
gene_sets <- split(fges$Gene, fges$Gene.signature)
bulkdata <- cpm(as.matrix(bulkdata), log = T)
bulkdata <- ExpressionSet(assayData = bulkdata)

bulkPar <- gsvaParam(exprData = bulkdata, geneSets = gene_sets, maxDiff=FALSE, kcdf = "Gaussian")
bulk_es <- gsva(bulkPar)

#Plot Heatmap
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
hmcol <- hmcol[length(hmcol):1]
pheatmap(exprs(bulk_es), col = hmcol, scale = "row",
         cluster_cols = T, cluster_rows = T,
         clustering_method = "ward.D2", angle_col = "45", fontsize_number = 12,
         fontsize = 12, height = 10)
