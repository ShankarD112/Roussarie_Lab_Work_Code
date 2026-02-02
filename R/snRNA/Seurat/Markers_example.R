

suppressPackageStartupMessages({
  library(Seurat)
  library(openxlsx)
})

############################################################
# PATHS (edit)
############################################################
sns_rds_path <- "path"
out_xlsx     <- "path"

############################################################
# Load Seurat object + quick UMAP plots
############################################################
SNSM <- readRDS(sns_rds_path)

DimPlot(SNSM, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(SNSM, reduction = "umap", group.by = "orig.ident", label = TRUE)

############################################################
# Define / load your subset
############################################################
# If SNSM_Exc_subset already exists in your environment, you can remove this block.
# Otherwise, make a subset here. Example:
# SNSM_Exc_subset <- subset(SNSM, subset = celltype == "Excitatory")

# Safety check
if (!exists("SNSM_Exc_subset")) {
  stop("SNSM_Exc_subset does not exist. Create it (subset from SNSM) before running markers.")
}

############################################################
# Prep SCT markers + FindAllMarkers
############################################################
SNSM_Exc_subset <- PrepSCTFindMarkers(SNSM_Exc_subset)

markers <- FindAllMarkers(
  object         = SNSM_Exc_subset,
  min.pct        = 0.25,
  logfc.threshold = 0
)

markers <- markers[markers$p_val_adj < 0.05, ]

############################################################
# Split by cluster + export to Excel (one sheet per cluster)
############################################################
marker_list <- split(markers, markers$cluster)
marker_list$All_Clusters <- markers

dir.create(dirname(out_xlsx), showWarnings = FALSE, recursive = TRUE)
write.xlsx(marker_list, file = out_xlsx, overwrite = TRUE)

cat("Wrote markers to:", out_xlsx, "\n")
