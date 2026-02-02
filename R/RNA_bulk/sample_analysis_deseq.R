# loading in the required packages
(library(DESeq2))
(library(tidyverse))
(library(readr))
(library(ggplot2))
(library(fgsea))
library(stats)
library(dplyr)


# Loading in files.
counts <- as.matrix(read.csv("path", row.name = "gene"))
id2gene <- read_delim("path", col_names = c("geneids", "genenames"))

print(head(counts))


coldata <- data.frame(samples = colnames(counts))

# Generic grouping: rename "groupA" / "groupB" to whatever matches your sample naming
coldata$Treatment <- ifelse(grepl("groupA", coldata$samples), "GroupA", "GroupB")
coldata$Treatment <- factor(coldata$Treatment, levels = c("GroupB", "GroupA"))

rownames(coldata) <- coldata$samples
coldata



# Example gene
plotCounts(dds, gene = "GENE_ID_1", intgroup = "Treatment")
# Example gene
plotCounts(dds, gene = "GENE_ID_2", intgroup = "Treatment")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~Treatment)
dds <- DESeq(dds)
res <- results(dds, contrast = c("Treatment", "GroupA", "GroupB"))

results <- res %>%
  as_tibble(rownames = "geneids") %>%
  left_join(id2gene, by = "geneids") %>%
  arrange(padj) %>%
  select(geneids, genenames, padj, pvalue, log2FoldChange)

print(results)
# Save results to CSV file
write.csv(results, "DESeq2_results.csv", row.names = FALSE)
filtered_results <- subset(results, padj <= 0.05)
filtered_results
summary(filtered_results$log2FoldChange)


resOrdered <- res[order(res$pvalue),] %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column()

resOrdered

library(ggplot2)
library(ggrepel)

vst_data <- varianceStabilizingTransformation(dds, blind = FALSE)

# Calculate PCA
pca_result <- prcomp(t(assay(vst_data)))

# Extract scores (samples' coordinates)
scores <- as.data.frame(pca_result$x)

# Create a data frame for plotting, including the treatment information
scores$sample <- rownames(scores)
scores$treatment <- coldata$Treatment[match(scores$sample, rownames(coldata))]

# Plot PCA biplot
ggplot(scores, aes(x = PC1, y = PC2, color = treatment, label = sample)) +
  geom_point(alpha = 0.8, size = 3) +
  geom_text_repel(aes(label = ifelse(abs(PC1) > 0.1 & abs(PC2) > 0.1, as.character(sample), "")),
                  size = 3, box.padding = 0.35, point.padding = 0.3) +
  xlab(paste("PC1: ", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance")) +
  ylab(paste("PC2: ", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance")) +
  ggtitle("PCA Biplot of RNA-seq Samples") +
  theme_minimal()



plot_log2fc <- results[results$padj < 0.05, ]

log2fc_plot <- ggplot(data = plot_log2fc, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.1, fill = "coral", color = "black") +
  labs(x = "Log 2 Fold Change", "Frequency",
       title = "Histogram of log2FoldChange Distributions")

print(log2fc_plot)


library(ggrepel)

# Define color based on the specified conditions
results$color <- with(results,
                      ifelse(padj < 0.05 & log2FoldChange >= 1, "Up",
                             ifelse(padj < 0.05 & log2FoldChange < -1, "Down",
                                    ifelse(padj < 0.05, "Significant", "Not Significant"))))

# Select top 10 most significant genes for labeling
top_genes <- results %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  slice_head(n = 10)

# Create the volcano plot with custom colors
volcano_plot <- ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Significant" = "black", "Down" = "blue", "Up" = "red", "Not Significant" = "gray")) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-value", title = "Volcano Plot of DE Genes") +
  theme_minimal() +
  geom_text_repel(data = top_genes, aes(label = genenames), box.padding = 0.35, point.padding = 0.2,
                  size = 3, max.overlaps = 10)

print(volcano_plot)

upregulated_genes <- results[results$color == "Up", ]
downregulated_genes <- results[results$color == "Down", ]

upregulated_genes <- upregulated_genes %>% select(geneids, genenames, padj, log2FoldChange)
downregulated_genes <- downregulated_genes %>% select(geneids, genenames, padj, log2FoldChange)

# Load gene sets using fgsea package
gene_sets <- fgsea::gmtPathways("path")

# Filter results for significantly differentially expressed genes
filtered_results <- subset(results, padj <= 0.05)

# Arrange results by log2FoldChange for use in GSEA analysis
rnk_list <- results %>%
  as_tibble() %>%
  arrange(desc(log2FoldChange)) %>%
  pull(log2FoldChange, name = genenames)

# Perform GSEA analysis
gsea_res <- fgsea(gene_sets, rnk_list, 15, 500)

# Select top pathways with positive NES
top_pos_nes <- gsea_res %>%
  arrange(desc(NES)) %>%
  select(pathway, padj, NES) %>%
  head(20)

# Select top pathways with negative NES
top_neg_nes <- gsea_res %>%
  arrange(NES) %>%
  select(pathway, padj, NES) %>%
  head(20)

print(top_pos_nes)
print(top_neg_nes)





