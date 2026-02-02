

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# Optional: source helper functions (comment out for portability)
# source("path")

############################################################
# PATHS (edit)
############################################################
base_dir  <- "path"   # directory containing per-sample Cell Ranger outputs
out_rds   <- "path"   # output RDS path

############################################################
# Samples (generic)
############################################################
samples <- paste0("sample", 1:13)

############################################################
# Load each sample as a Seurat object
############################################################
seurat_objects <- vector("list", length(samples))
names(seurat_objects) <- samples

for (sample in samples) {
  data_dir <- file.path(base_dir, sample, "outs", "filtered_feature_bc_matrix")

  if (!dir.exists(data_dir)) {
    stop("Missing directory: ", data_dir)
  }

  counts <- Read10X(data_dir = data_dir)
  seurat_objects[[sample]] <- CreateSeuratObject(counts = counts, project = sample)
}

############################################################
# Merge into one Seurat object
############################################################
combined <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1]
)

############################################################
# QC: mitochondrial percent and filtering
############################################################
# NOTE: pattern depends on gene naming convention:
# - mouse often: "^mt-"
# - human often: "^MT-"
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

combined_filtered <- subset(combined, subset = percent.mt < 5)

############################################################
# SCTransform + downstream
############################################################
combined_filtered <- SCTransform(
  combined_filtered,
  vars.to.regress = "nFeature_RNA",
  verbose = FALSE,
  seed.use = 1448145
)

combined_filtered <- RunPCA(combined_filtered, npcs = 50, verbose = FALSE)
combined_filtered <- RunUMAP(combined_filtered, dims = 1:50, verbose = FALSE)
combined_filtered <- FindNeighbors(combined_filtered, dims = 1:50, verbose = FALSE)
combined_filtered <- FindClusters(combined_filtered, resolution = 1, verbose = FALSE)

############################################################
# Save
############################################################
dir.create(dirname(out_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(combined_filtered, file = out_rds)

cat("Saved filtered + SCT Seurat object to:", out_rds, "\n")
