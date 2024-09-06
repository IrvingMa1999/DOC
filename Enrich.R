################## Generate Background File ###########################
library(KEGGREST)

# Read in functional prediction data
KEGG <- read.csv("functional_prediction.txt", row.names = 1, header = TRUE, sep = "\t")

# Subset the data
KEGG_4612_end <- KEGG[5475:9166, 1:1013]

# Extract enzyme IDs
enzyme_ids <- rownames(KEGG_4612_end)

# Initialize dataframe to store enzyme and pathway associations
all_enzyme_pathways <- data.frame(enzyme_id = character(), pathway_id = character(), stringsAsFactors = FALSE)

# Iterate over enzyme IDs and query associated pathways
for (enzyme_id in enzyme_ids) {
  pathways <- keggLink("pathway", enzyme_id)
  df <- data.frame(
    enzyme_id = rep(enzyme_id, length(pathways)),
    pathway_id = pathways,
    stringsAsFactors = FALSE
  )
  all_enzyme_pathways <- rbind(all_enzyme_pathways, df)
  Sys.sleep(1)  # Pause to avoid overloading the KEGG server
}

# Save the enzyme-pathway associations to a CSV file
write.csv(all_enzyme_pathways, "last_pathway.csv")

# Read in the total pathway gene mapping
gene_pathway <- read.csv("total_pathway.csv", header = TRUE)

############## Enrichment Analysis ########################

# Set working directory
setwd("W:/16S_NEW/enrich/subgroup/")

# Read in differential count data and sample information
diffcount <- read.csv("Tropical_Climate_function.csv", header = TRUE, row.names = 1)
sampleTable <- read.csv("Tropical_Climate_sampel_function.csv", header = TRUE, row.names = 1)

# Remove rows with NA values
diffcount <- na.omit(diffcount)

# Convert differential count data to a matrix and multiply by 10,000 for scaling
diffcount <- as.matrix(diffcount) * 10000

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(diffcount), colData = sampleTable, design = ~group)

# Normalize the DESeq2 dataset
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "High", "Low"))

# Convert results to a dataframe
diff_res <- as.data.frame(res)

# Filter results by p-value
diff_p <- diff_res[order(diff_res$pvalue),]
diff_gene_deseq2_p <- subset(diff_p, pvalue < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))

# Save significant differential expression results to a CSV file
write.csv(diff_gene_deseq2_p, file = "deseq_tropical.csv")

# Extract background and significant genes
background_genes <- rownames(diffcount)
significant_genes <- rownames(diff_gene_deseq2_p)

# Perform enrichment analysis
library(clusterProfiler)
gene2kegg <- gene_pathway
enrich_result <- enricher(significant_genes,
                          TERM2GENE = gene2kegg[, c("KEGG_Pathway", "query_name")],
                          TERM2NAME = NULL,  # Set to NULL if pathway names are not available
                          universe = background_genes,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

# Extract and save enrichment results
enrich_res <- enrich_result@result
write.csv(enrich_res, "enrichkegg_res.csv")