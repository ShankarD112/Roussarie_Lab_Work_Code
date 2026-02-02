# Reusable single-cell helper functions (Seurat + friends).
# Usage: source("/path/single_cell.R") in your script.


# Load necessary library


suppressMessages({
  suppressWarnings({
    library(sctransform)
    library(Seurat)
    library(glmGamPoi)
    library(Matrix)
    library(SeuratDisk)
    library(sctransform)
    library(glmGamPoi)
    library(dplyr)
    library(ggplot2)
    library(reticulate)
    library(DoubletFinder)
    library(purrr)
    library(ggridges)
    library(writexl)
    library(stringr)
    library(tibble)
    library(readxl)
    library(readr)
    library(tidyr)
    library(pheatmap)
    library(knitr)
    library(openxlsx)
    library(scDblFinder)
    library(data.table)
    library(tximport) 
    library(scuttle)
    library(SingleCellExperiment)
  })
})


# Define a reusable function to load the samples and calculate the mito percentage
load_samples <- function(base_dir, sample_names, pattern = "^MT-", min_cells = 3, min_features = 200) {
  # Load required library
  library(Seurat)
  
  # Initialize an empty list to store the Seurat objects
  sample_objects <- list()
  
  # Loop through sample names
  for (s in sample_names) {
    # Path to the filtered feature matrix
    data_dir <- file.path(base_dir, s, "outs", "filtered_feature_bc_matrix")
    
    # Check if the directory exists
    if (!dir.exists(data_dir)) {
      warning(paste("Directory not found for sample:", s, "Skipping..."))
      next
    }
    
    # Load the data
    counts <- Read10X(data.dir = data_dir)
    
    # Append sample name to barcodes
    colnames(counts) <- paste0(colnames(counts), "_", s)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = counts, project = s, min.cells = min_cells, min.features = min_features)
    
    # Calculate mitochondrial percentage
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = pattern)
    
    # Add the Seurat object to the list
    sample_objects[[s]] <- seurat_obj
  }
  
  # Return the list of Seurat objects
  return(sample_objects)
}



# Define a reusable plotting function
generate_seurat_plots <- function(seurat_object, reduction = "umap", group_by_list = NULL, feature_list = NULL, label = TRUE) {
  # Check for group_by_list and plot DimPlots
  if (!is.null(group_by_list)) {
    for (group_by in group_by_list) {
      print(DimPlot(seurat_object, reduction = reduction, label = label, group.by = group_by,
                    na.value = "grey90"))
    }
  }
  
  # Check for feature_list and plot FeaturePlots
  if (!is.null(feature_list)) {
    for (feature in feature_list) {
      print(FeaturePlot(seurat_object, features = feature))
    }
  }
}


# Define a reusable function to save Seurat plots
save_seurat_plots <- function(seurat_object, output_dir, reduction = "umap", group_by_list = NULL, feature_list = NULL, label = TRUE, plot_width = 8, plot_height = 6, dpi = 500) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save DimPlots
  if (!is.null(group_by_list)) {
    for (group_by in group_by_list) {
      plot <- DimPlot(seurat_object, reduction = reduction, label = label, group.by = group_by)
      filename <- file.path(output_dir, paste0("UMAP_", group_by, ".png"))
      ggsave(filename = filename, plot = plot, width = plot_width, height = plot_height, dpi = dpi)
    }
  }
  
  # Save FeaturePlots
  if (!is.null(feature_list)) {
    for (feature in feature_list) {
      plot <- FeaturePlot(seurat_object, features = feature)
      filename <- file.path(output_dir, paste0("FeaturePlot_", feature, ".png"))
      ggsave(filename = filename, plot = plot, width = plot_width, height = plot_height, dpi = dpi)
    }
  }
}


# Define a reusable function to save cluster markers
save_markers_to_excel <- function(markers, output_path) {
  
  # Split markers by cluster
  marker_list <- split(markers, markers$cluster)
  
  # Add a combined sheet to the list
  marker_list$All_Clusters <- markers
  
  # Write markers to an Excel file with multiple sheets
  write.xlsx(marker_list, file = output_path)
  
  message("Markers saved to Excel file at: ", output_path)
}

# Define a reusable function to visualize genes

PlotGeneList <- function(
    seurat_obj,
    gene_list,
    assay       = "RNA",
    reduction   = "umap",
    split.by = NULL
) {
  # Set the default assay
  DefaultAssay(seurat_obj) <- assay
  
  # Loop over each gene and print a FeaturePlot
  for (gene in gene_list) {
    print(
      FeaturePlot(
        object     = seurat_obj,
        features   = gene,
        reduction  = reduction,
        split.by = split.by,
        cols       = c("yellow", "blue")  # color scale: grey -> red
      ) #+
        #ggtitle(gene)  # label the plot with the gene name
    )
  }
}

# Function to map mouse genes to human genes and remove unmapped genes
load_orthologs_names <- function(
    sample_names, 
    dir_path, 
    orthologs_path,
    include_unmapped = TRUE
) {
  # Load the orthologs table
  orthologs <- read.csv(orthologs_path)
  orthologs <- orthologs[!duplicated(orthologs$mouse_gene), ]
  orthologs <- orthologs[complete.cases(orthologs), ]
  
  # Map mouse genes to human genes
  gene_mapping <- setNames(orthologs$human_gene, orthologs$mouse_gene)
  
  # Initialize a list to store Seurat objects
  seurat_objects <- list()
  
  for (sample_name in sample_names) {
    # Construct the file path
    file_path <- file.path(dir_path, sample_name, "outs/filtered_feature_bc_matrix.h5")
    
    # Read the matrix data
    matrix_data <- Read10X_h5(file_path)
    
    # Map and update row names
    new_row_names <- gene_mapping[rownames(matrix_data)]
    
    if (include_unmapped) {
      # Replace row names where mappings exist, leave others unchanged
      rownames(matrix_data)[!is.na(new_row_names)] <- new_row_names[!is.na(new_row_names)]
    } else {
      # Filter to only include genes with valid mappings
      matrix_data <- matrix_data[!is.na(new_row_names), ]
      rownames(matrix_data) <- new_row_names[!is.na(new_row_names)]
    }
    
    
    # Create a Seurat object for the sample
    seurat_object <- CreateSeuratObject(counts = matrix_data, project = sample_name)
    seurat_objects[[sample_name]] <- seurat_object
  }
  
  # Merge all Seurat objects into one
  merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_objects)
  
  return(merged_seurat)
}


# Function to check the presence of a gene in the seurat object
check_gene_presence <- function(seurat_obj, gene) {
  if (!(gene %in% rownames(seurat_obj))) {
    message(paste("❌ Gene", gene, "is NOT found in the Seurat object."))
    return(NULL)
  }
  
  # Extract gene expression values
  gene_counts <- FetchData(seurat_obj, vars = gene)
  
  # Count number of cells where the gene is expressed (nonzero values)
  num_cells_expressed <- sum(gene_counts > 0)
  
  message(paste("✅ Gene", gene, "is present in the Seurat object."))
  message(paste("The gene is expressed in", num_cells_expressed, "cells out of", ncol(seurat_obj), "total cells."))
  
  return(num_cells_expressed)
}

# Define a reusable function to visualize all genes

PlotGeneLista <- function(
    seurat_obj,
    gene_list,
    assay       = "SCT",
    reduction   = "umap"
) {
  # Set the default assay
  DefaultAssay(seurat_obj) <- assay
  
  # Extract expression matrix for the given genes
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = "data")[gene_list, ]
  
  # Ensure all genes are present in the dataset
  missing_genes <- setdiff(gene_list, rownames(expr_matrix))
  if (length(missing_genes) > 0) {
    warning(paste("The following genes are not found in the dataset and will be ignored:", paste(missing_genes, collapse = ", ")))
  }
  
  # Identify cells where ALL genes in the list are expressed (nonzero)
  cells_to_color <- colSums(expr_matrix > 0) == length(gene_list)
  
  # Create a new metadata column to store binary classification
  seurat_obj$GeneList_Expression <- ifelse(cells_to_color, "Expressed", "Not Expressed")
  
  # Plot UMAP with binary coloring
  print(
    DimPlot(
      seurat_obj,
      group.by   = "GeneList_Expression",
      reduction  = reduction#,
      #cols       = c("blue", "yellow")  # Not expressed -> grey, Expressed -> red
    ) +
      ggtitle(paste("Cells expressing all genes in list"))
  )
}




PlotGeneList_CoGE <- function(
    seurat_obj,
    gene_list,
    assay = "SCT",
    reduction = "umap"
) {
  DefaultAssay(seurat_obj) <- assay
  
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = "data")[gene_list, , drop = FALSE]
  
  cell_labels <- rep("None", ncol(expr_matrix))
  
  for (gene in gene_list) {
    expressed <- expr_matrix[gene, ] > 0
    cell_labels[expressed & cell_labels == "None"] <- gene
  }
  
  seurat_obj$GeneExpressionLabel <- factor(cell_labels, levels = c(gene_list, "None"))
  
  p <- DimPlot(
    seurat_obj,
    group.by = "GeneExpressionLabel",
    reduction = reduction,
    cols = c(scales::hue_pal()(length(gene_list)), "lightgrey")  # grey for "None"
  ) +
    ggtitle("Cells colored by gene expression (first match)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}
