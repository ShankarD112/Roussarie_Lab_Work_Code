#!/usr/bin/env Rscript

# Patch-seq WGCNA workflow (cleaned, GitHub-ready)
# - Removed Rmd chunk markers
# - Replaced HPC paths with "path"
# - Kept analysis flow intact
# - NOTE: This script assumes Seurat/tidyverse/readxl/etc. are installed
# - Edit the PATHS section before running

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(openxlsx)
  library(DESeq2)
  library(WGCNA)
  library(flashClust)
  library(curl)
  library(dendextend)
})

options(stringsAsFactors = FALSE)
WGCNA::allowWGCNAThreads()

############################################################
# PATHS (edit)
############################################################

single_cell_functions <- "path"

tpm_csv_path          <- "path"
seurat_rds_path       <- "path"
metadata_month_path   <- "path"
metadata_path         <- "path"
mapping_path          <- "path"

out_dir_wgcna         <- "path"
dir.create(out_dir_wgcna, showWarnings = FALSE, recursive = TRUE)

############################################################
# 0. Source helper functions (optional)
############################################################

# If you don't want to depend on an external file for GitHub,
# either vendor/copy the functions into this repo or comment this out.
if (file.exists(single_cell_functions)) {
  source(single_cell_functions)
}

############################################################
## 1. Read TPM matrix (introns only, no mt genes)
############################################################

intron_only_no_mt <- read.csv(
  tpm_csv_path,
  row.names   = 1,
  check.names = FALSE
)

############################################################
## 2. Read Seurat object + basic metadata
############################################################

cogent_merged <- readRDS(seurat_rds_path)
DefaultAssay(cogent_merged) <- "RNA"

metadata_selected <- cogent_merged@meta.data %>%
  dplyr::select(Sample, Intron_to_Unique_ratio, percent.mt)

print(head(metadata_selected))

############################################################
## 3. Read external metadata (month, main ephys, mapping)
############################################################

metadata_month <- read_excel(metadata_month_path, sheet = 1)

metadata <- read_excel(metadata_path, sheet = 1) %>%
  mutate(Sample = paste0("AM", Sample)) %>%
  filter(Sample != "AMNA")

exclude_cols <- c(
  "Sample", "Location", "MECvsLEC", "Final_Determination", "Morphology",
  "Do_we_see_cell_in_L2EC", "Grik1", "Gpc5", "Tenm3",
  "Genotype", "Sex", "Age_days"
)

metadata <- metadata %>%
  mutate(
    across(
      .cols = !all_of(exclude_cols),
      .fns = ~ {
        if (is.character(.) || is.factor(.)) {
          suppressWarnings(as.numeric(as.character(.)))
        } else {
          .
        }
      }
    )
  )

mapping <- read_excel(mapping_path, sheet = 1)

samples_metadata <- metadata$Sample
samples_mapping  <- mapping$Sample

common_samples <- intersect(samples_metadata, samples_mapping)
cat("Common samples:", length(common_samples), "\n")

metadata_common <- metadata[metadata$Sample %in% common_samples, ]
mapping_common  <- mapping[mapping$Sample %in% common_samples, ]

metadata_common <- metadata_common[order(metadata_common$Sample), ]
mapping_common  <- mapping_common[order(mapping_common$Sample), ]

metadata_common <- inner_join(metadata_common, mapping_common, by = "Sample")

metadata_common <- metadata_common %>%
  mutate(
    match_mismatch = case_when(
      MMC_Region == Final_Determination ~ "all_match",
      MMC_Region != Final_Determination ~ "mismatch_MMC_Region"
    )
  )

metadata_common <- inner_join(metadata_common, metadata_month, by = "Sample")
metadata_common$Final_Determination <- toupper(metadata_common$Final_Determination)

print(metadata_common)

############################################################
## 4. Restrict to LEC cells & align with TPM matrix
############################################################

lec_meta <- metadata_common %>%
  filter(
    Final_Determination == "LEC",
    Genotype %in% c("APP", "WT")
  )

lec_cells <- intersect(lec_meta$Sample, colnames(intron_only_no_mt))
print(lec_cells)

lec_counts <- intron_only_no_mt[, lec_cells, drop = FALSE]

meta_lec <- metadata_common %>%
  filter(
    Final_Determination == "LEC",
    Sample %in% colnames(lec_counts)
  ) %>%
  dplyr::filter(Genotype %in% c("WT", "APP"))

keep_cells <- intersect(meta_lec$Sample, colnames(lec_counts))
cat("WT+APP cells retained:", length(keep_cells), "\n")

lec_counts <- lec_counts[, keep_cells, drop = FALSE]

print(table(meta_lec$Genotype))

############################################################
## 5. Filter low expressed genes
############################################################

n_samples <- ncol(lec_counts)
keep_genes <- rowSums(lec_counts > 0) >= (n_samples * 0.25)

cat("Genes retained:", sum(keep_genes), "out of", nrow(intron_only_no_mt), "\n")

lec_counts_filt <- lec_counts[keep_genes, ]
print(dim(lec_counts_filt))

############################################################
## 6. Variance stabilizing transform (DESeq2)
############################################################

sample_info <- data.frame(row.names = colnames(lec_counts_filt))

dds <- DESeqDataSetFromMatrix(
  countData = lec_counts_filt,
  colData   = sample_info,
  design    = ~ 1
)

dds <- estimateSizeFactors(dds)
vst_dds <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_dds)

# WGCNA expects samples x genes
vst_mat_df <- as.data.frame(t(vst_mat))

vst_mat_df_gsg <- goodSamplesGenes(vst_mat_df)
cat("goodSamplesGenes allOK:", vst_mat_df_gsg$allOK, "\n")

############################################################
## 7. Sample clustering + trait color bars
############################################################

sampleTree <- hclust(dist(vst_mat_df), method = "average")

par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")

dend_no_reorder <- as.dendrogram(sampleTree)
dend_no_reorder <- rotate(dend_no_reorder, order = 1:length(sampleTree$labels))

par(cex = 0.6, mar = c(10, 4, 2, 0))
plot(dend_no_reorder, horiz = TRUE, main = "Sample clustering (no leaf reordering)")

meta <- cogent_merged@meta.data

percent_mt_by_sample <- meta %>%
  group_by(Sample) %>%
  summarise(percent.mt = mean(percent.mt, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

percent_mt <- percent_mt_by_sample$percent.mt
names(percent_mt) <- percent_mt_by_sample$Sample

percent_mt_ordered <- percent_mt[sampleTree$labels]

sampleColors <- numbers2colors(percent_mt_ordered, signed = FALSE, naColor = "grey")
par(cex = 0.6, mar = c(0, 4, 2, 0))
plotDendroAndColors(sampleTree, sampleColors, groupLabels = "percent.mt",
                    main = "Sample clustering with percent.mt")

intron_ratio_by_sample <- meta %>%
  group_by(Sample) %>%
  summarise(Intron_to_Unique_ratio = mean(Intron_to_Unique_ratio, na.rm = TRUE), .groups = "drop") %>%
  as.data.frame()

intron_ratio <- intron_ratio_by_sample$Intron_to_Unique_ratio
names(intron_ratio) <- intron_ratio_by_sample$Sample

intron_ratio_ordered <- intron_ratio[sampleTree$labels]
intronColors <- numbers2colors(intron_ratio_ordered, signed = FALSE, naColor = "grey")

par(cex = 0.6, mar = c(1, 8, 3, 1))
plotDendroAndColors(sampleTree, intronColors,
                    groupLabels = "Intron_to_Unique_ratio",
                    main = "Sample clustering with Intron_to_Unique_ratio",
                    cex.colorLabels = 0.4)

############################################################
## 8. Pick soft threshold
############################################################

spt <- pickSoftThreshold(vst_mat_df)

par(mar = c(5, 5, 4, 2))
plot(spt$fitIndices[, 1], spt$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(spt$fitIndices[, 1], spt$fitIndices[, 2], labels = spt$fitIndices[, 1], col = "red")
abline(h = 0.80, col = "red")

plot(spt$fitIndices[, 1], spt$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(spt$fitIndices[, 1], spt$fitIndices[, 5], labels = spt$fitIndices[, 1], col = "red")

softPower <- 3  # chosen power

############################################################
## 9. Network construction + modules
############################################################

adjacency <- adjacency(vst_mat_df, power = softPower)
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM

geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

sizeGrWindow(12, 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

Modules <- cutreeDynamic(
  dendro = geneTree,
  distM = TOM.dissimilarity,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)

ModuleColors <- labels2colors(Modules)

plotDendroAndColors(geneTree, ModuleColors, "Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

MElist <- moduleEigengenes(vst_mat_df, colors = ModuleColors)
MEs <- MElist$eigengenes

ME.dissimilarity <- 1 - cor(MEs, use = "complete")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")

par(mar = c(0, 4, 2, 0), cex = 0.6)
plot(METree)
abline(h = .25, col = "red")

merge <- mergeCloseModules(vst_mat_df, ModuleColors, cutHeight = .25)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors),
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors (original vs merged)")

ME.dissimilarity.post <- 1 - cor(mergedMEs, use = "complete")
METree.post <- hclust(as.dist(ME.dissimilarity.post), method = "average")

par(mar = c(0, 4, 2, 0), cex = 0.6)
plot(METree.post, main = "Merged module eigengene clustering", sub = "", xlab = "")
abline(h = 0.25, col = "red")

############################################################
## 10. Build datTraits from meta_lec
############################################################

traitData <- meta_lec %>%
  dplyr::filter(Sample %in% rownames(vst_mat_df))

traitData <- traitData[match(rownames(vst_mat_df), traitData$Sample), ]

numeric_cols <- names(traitData)[sapply(traitData, is.numeric)]
traitData$Genotype_APP <- ifelse(traitData$Genotype == "APP", 1, 0)
numeric_cols <- c(numeric_cols, "Genotype_APP")

datTraits <- traitData[, c("Sample", numeric_cols)]
rownames(datTraits) <- datTraits$Sample

stopifnot(all(rownames(datTraits) == rownames(vst_mat_df)))

############################################################
## 11. Module–trait correlations + heatmap
############################################################

nSamples <- nrow(vst_mat_df)

module.trait.correlation <- cor(mergedMEs, datTraits, use = "pairwise.complete.obs")
module.trait.Pvalue <- corPvalueStudent(module.trait.correlation, nSamples)

textMatrix <- paste(
  signif(module.trait.correlation, 2), "\n(",
  signif(module.trait.Pvalue, 1), ")", sep = ""
)
dim(textMatrix) <- dim(module.trait.correlation)

par(mar = c(6, 8.5, 3, 1))
labeledHeatmap(
  Matrix = module.trait.correlation,
  xLabels = names(datTraits),
  yLabels = names(mergedMEs),
  ySymbols = names(mergedMEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.4,
  zlim = c(-1, 1),
  main = "Module–trait relationships"
)

############################################################
## 12. Save module–trait r/p side-by-side to Excel
############################################################

traits <- colnames(module.trait.correlation)
combined_df <- data.frame(row.names = rownames(module.trait.correlation))

for (tr in traits) {
  combined_df[[paste0(tr, "_r")]] <- module.trait.correlation[, tr]
  combined_df[[paste0(tr, "_p")]] <- module.trait.Pvalue[, tr]
}

out_xlsx_module_trait <- file.path(out_dir_wgcna, "Module_Trait_Relationships.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "r_p_side_by_side")
writeData(wb, sheet = "r_p_side_by_side", x = combined_df, rowNames = TRUE)
saveWorkbook(wb, out_xlsx_module_trait, overwrite = TRUE)

############################################################
## 13. Save large PDF heatmap
############################################################

pdf_out1 <- file.path(out_dir_wgcna, "Patch_seq_WGCNA_LEC_Module_trait_relationships_big.pdf")

pdf(file = pdf_out1, width = 40, height = 32, pointsize = 18)
par(mar = c(16, 28, 6, 10))
labeledHeatmap(
  Matrix = module.trait.correlation,
  xLabels = names(datTraits),
  yLabels = names(mergedMEs),
  ySymbols = names(mergedMEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.25,
  cex.lab = 0.7,
  xLabelsAngle = 90,
  zlim = c(-1, 1),
  main = "Module–trait relationships"
)
dev.off()

############################################################
## 14. Optional: filter out certain frequency/time traits and re-plot/save
############################################################

all_numeric_meta <- colnames(datTraits)
freq_regex <- "(80|130|180|230|280|330)"

to_drop <- all_numeric_meta[
  grepl(paste0("^Low_", freq_regex, "_freq_(1st|2nd|3rd)$"), all_numeric_meta) |
    grepl(paste0("^", freq_regex, "_freq_(1st|2nd)_half$"), all_numeric_meta) |
    grepl("^Adapt_ratio_3rd_1st_av\\.freq", all_numeric_meta) |
    grepl(paste0("^", freq_regex, "_amp_1st(\\.|$)"), all_numeric_meta) |
    grepl(paste0("^", freq_regex, "_amp_last(\\.|$)"), all_numeric_meta) |
    grepl("^time1st(\\.|$)", all_numeric_meta) |
    grepl("^timelast(\\.|$)", all_numeric_meta)
]

keep_traits <- setdiff(all_numeric_meta, to_drop)
datTraits_keep <- datTraits[, keep_traits, drop = FALSE]

module.trait.correlation_f <- cor(mergedMEs, datTraits_keep, use = "pairwise.complete.obs")
module.trait.Pvalue_f <- corPvalueStudent(module.trait.correlation_f, nSamples)

textMatrix_f <- paste(
  signif(module.trait.correlation_f, 2), "\n(",
  signif(module.trait.Pvalue_f, 1), ")", sep = ""
)
dim(textMatrix_f) <- dim(module.trait.correlation_f)

pdf_out2 <- file.path(out_dir_wgcna, "Patch_seq_WGCNA_LEC_Module_trait_relationships_big_FILTERED.pdf")

pdf(file = pdf_out2, width = 34, height = 28, pointsize = 18)
par(mar = c(14, 26, 6, 8))
labeledHeatmap(
  Matrix = module.trait.correlation_f,
  xLabels = colnames(datTraits_keep),
  yLabels = names(mergedMEs),
  ySymbols = names(mergedMEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix_f,
  setStdMargins = FALSE,
  cex.text = 0.3,
  cex.lab = 0.7,
  xLabelsAngle = 90,
  zlim = c(-1, 1),
  main = "Module–trait relationships (freq vars removed)"
)
dev.off()

############################################################
## 15. Save module gene lists to Excel
############################################################

genes <- colnames(vst_mat_df)

gene_module_df <- data.frame(
  Gene = genes,
  Module = ModuleColors,
  stringsAsFactors = FALSE
)

xlsx_gene_lists <- file.path(out_dir_wgcna, "WGCNA_Module_Gene_Lists.xlsx")

wb2 <- createWorkbook()
addWorksheet(wb2, "All_Genes")
writeData(wb2, "All_Genes", gene_module_df)

module_colors_unique <- sort(unique(ModuleColors))
for (mod in module_colors_unique) {
  mod_genes <- gene_module_df$Gene[gene_module_df$Module == mod]
  df_mod <- data.frame(Gene = mod_genes, stringsAsFactors = FALSE)
  sheet_name <- substr(paste0("Module_", mod), 1, 31)
  addWorksheet(wb2, sheet_name)
  writeData(wb2, sheet_name, df_mod)
}
saveWorkbook(wb2, xlsx_gene_lists, overwrite = TRUE)

############################################################
## 16. Plot MM vs GS helper + example usage
############################################################

plotMMvsGS <- function(module,
                       traitName,
                       datExpr,
                       MEs,
                       colors,
                       datTraits,
                       useAbs = TRUE,
                       ...) {
  library(WGCNA)

  commonSamples <- intersect(rownames(datExpr), rownames(datTraits))
  datExpr   <- datExpr[commonSamples, , drop = FALSE]
  MEs       <- MEs[commonSamples, , drop = FALSE]
  datTraits <- datTraits[commonSamples, , drop = FALSE]

  if (!traitName %in% colnames(datTraits)) {
    stop("Trait '", traitName, "' not found in datTraits.")
  }

  traitDF <- as.data.frame(datTraits[[traitName]])
  rownames(traitDF) <- commonSamples
  names(traitDF) <- traitName

  modNames <- substring(names(MEs), 3)
  nSamples <- nrow(datExpr)

  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)

  geneTraitSignificance <- as.data.frame(cor(datExpr, traitDF, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) <- paste0("GS.", traitName)
  names(GSPvalue) <- paste0("p.GS.", traitName)

  colIdx <- match(module, modNames)
  if (is.na(colIdx)) stop("Module '", module, "' not found in modNames.")

  moduleGenes <- colors == module

  x <- as.numeric(geneModuleMembership[moduleGenes, colIdx])
  y <- as.numeric(geneTraitSignificance[moduleGenes, 1])

  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]

  if (useAbs) {
    x <- abs(x); y <- abs(y)
  }

  par(mar = c(5, 5, 4, 2))
  verboseScatterplot(
    x, y,
    xlab = paste("Module membership in", module, "module"),
    ylab = paste("Gene significance for", traitName),
    main = "Module membership vs. gene significance",
    col = module,
    ...
  )

  fit <- lm(y ~ x)
  xs <- seq(min(x), max(x), length.out = 200)
  pred <- predict(fit, newdata = data.frame(x = xs), interval = "confidence")

  polygon(
    c(xs, rev(xs)),
    c(pred[, "lwr"], rev(pred[, "upr"])),
    border = NA,
    col = grDevices::adjustcolor("grey50", alpha.f = 0.3)
  )
  lines(xs, pred[, "fit"], lwd = 2)

  r <- cor(x, y, use = "p")
  p <- corPvalueStudent(r, length(x))

  invisible(list(
    MM = geneModuleMembership,
    MM.p = MMPvalue,
    GS = geneTraitSignificance,
    GS.p = GSPvalue,
    moduleGenes = moduleGenes,
    cor = r,
    p = p
  ))
}

# Example calls (edit module/trait names as needed)
# plotMMvsGS("tan", "Sag_mV_Manually_measured_", vst_mat_df, mergedMEs, mergedColors, datTraits)

############################################################
## 17. Save core objects to an RData file
############################################################

save_file <- file.path(out_dir_wgcna, "Patch_seq_WGCNA_LEC_core.RData")

Genotype_APP <- as.data.frame(datTraits$Genotype_APP)
names(Genotype_APP) <- "Genotype_APP"

# Compute membership/significance once for saving (for downstream use)
modNames <- substring(names(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(vst_mat_df, mergedMEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(vst_mat_df, Genotype_APP, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(Genotype_APP), sep = "")
names(GSPvalue) <- paste("p.GS.", names(Genotype_APP), sep = "")

save(
  vst_mat_df,
  lec_counts_filt,
  keep_genes,
  meta_lec,
  traitData,
  datTraits,
  datTraits_keep,
  keep_traits,
  to_drop,
  softPower,
  spt,
  sampleTree,
  geneTree,
  Modules,
  ModuleColors,
  mergedColors,
  MElist,
  MEs,
  mergedMEs,
  module.trait.correlation,
  module.trait.Pvalue,
  gene_module_df,
  geneModuleMembership,
  MMPvalue,
  geneTraitSignificance,
  GSPvalue,
  adjacency,
  TOM,
  TOM.dissimilarity,
  Genotype_APP,
  file = save_file
)

cat("Saved core objects to:", save_file, "\n")

############################################################
## 18. Extra: eigengene network plots (optional)
############################################################

MET <- orderMEs(mergedMEs)

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0, 4, 2, 0),
                      plotHeatmaps = FALSE)

par(cex = 1.0, mar = c(1, 1, 1, 1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5, 5, 2, 2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

############################################################
## 19. Save Gene Module Membership tables (optional)
############################################################

MM_table <- cbind(geneModuleMembership, MMPvalue)

out_path_all <- file.path(out_dir_wgcna, "Gene_Module_Membership.xlsx")
out_dir_modules <- file.path(out_dir_wgcna, "Gene_module_membership_per_module")
dir.create(out_dir_modules, showWarnings = FALSE, recursive = TRUE)

wb3 <- createWorkbook()
addWorksheet(wb3, "ModuleMembership")
writeData(wb3, sheet = "ModuleMembership", x = MM_table, rowNames = TRUE)
saveWorkbook(wb3, out_path_all, overwrite = TRUE)

module_names <- gsub("^MM", "", grep("^MM", colnames(MM_table), value = TRUE))

for (mod in module_names) {
  mm_col <- paste0("MM", mod)
  p_col  <- paste0("p.MM", mod)
  module_df <- MM_table[, c(mm_col, p_col), drop = FALSE]

  out_file <- file.path(out_dir_modules, paste0("Gene_Module_Membership_", mod, ".xlsx"))

  wb_mod <- createWorkbook()
  sheet_name <- substr(mod, 1, 31)
  addWorksheet(wb_mod, sheet_name)
  writeData(wb_mod, sheet = sheet_name, x = module_df, rowNames = TRUE)
  saveWorkbook(wb_mod, out_file, overwrite = TRUE)
}

cat("Finished.\n")
