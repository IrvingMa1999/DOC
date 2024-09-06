# Load necessary libraries and install dependencies
library(DESeq2)
library(BiocManager)
BiocManager::install("apeglm")
library(apeglm)

# Read in the count matrix and sample information
count_matrix <- read.csv('otu.csv', header=TRUE, row.names=1)
sample_info <- read.csv('Sample_info.csv', header=TRUE, row.names=1)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Treatment)
dds

# Filter low expression genes
keep <- rowSums(counts(dds)) >= 6
dds <- dds[keep,]

# Perform differential expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast=c("Treatment", "A", "N"))
res
resultsNames(dds)

# Log fold change shrinkage
resLFC <- lfcShrink(dds, coef="A_vs_B", type="apeglm")
resLFC

# Order results by log2FoldChange
resOrdered <- res[order(res$log2FoldChange),]
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)

# Convert results to data frame and filter significant results
res <- as.data.frame(res)
resSig <- subset(res, padj < 0.1)
resSig

# Save significant results to CSV
write.csv(as.data.frame(resSig), file="results05.csv")