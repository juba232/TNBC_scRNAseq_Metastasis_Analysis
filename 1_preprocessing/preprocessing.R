# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


# Load the count matrix
counts <- read.csv("counts_rsem.csv.gz", row.names = 1, check.names = FALSE)

# Load metadata
meta <- read.csv("meta_data.csv", row.names = 1, check.names = FALSE)

# Create Seurat object with initial filtering
seurat <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  min.cells = 10,      # Keep genes expressed in ≥10 cells
  min.features = 200   # Keep cells with ≥200 genes
)

# Add basic QC metrics
seurat$log10GenesPerUMI <- log10(seurat$nFeature_RNA) / log10(seurat$nCount_RNA)
seurat$patient <- seurat$patient

# Calculate complexity (genes per UMI)
seurat$complexity <- seurat$nFeature_RNA / seurat$nCount_RNA

# Set QC thresholds based on percentiles
qc_metrics <- data.frame(
  nCount_RNA = seurat$nCount_RNA,
  nFeature_RNA = seurat$nFeature_RNA,
  complexity = seurat$complexity
)

# Set appropriate thresholds based on percentiles
corrected_count_threshold_low <- 3000    # ~5th percentile
corrected_count_threshold_high <- 2500000 # ~95th percentile  
corrected_feature_threshold_low <- 435   # ~5th percentile
corrected_complexity_threshold <- 0.0015 # ~5th percentile

# Apply corrected filtering
cells_keep_corrected <- seurat$nCount_RNA >= corrected_count_threshold_low & 
  seurat$nCount_RNA <= corrected_count_threshold_high &
  seurat$nFeature_RNA >= corrected_feature_threshold_low &
  seurat$complexity >= corrected_complexity_threshold

cat("Cells with corrected thresholds:", sum(cells_keep_corrected), "\n")

seurat_filtered <- subset(seurat, cells = colnames(seurat)[cells_keep_corrected])
cat("Successfully filtered to", ncol(seurat_filtered), "cells\n")

# Check what percentage we kept
cat("Percentage of cells retained:", round(ncol(seurat_filtered)/ncol(seurat)*100, 1), "%\n")

# Create QC plots with appropriate thresholds
library(ggplot2)
library(patchwork)

# Figure 1A: Count distribution by patient with corrected thresholds
p1 <- VlnPlot(seurat, features = "nCount_RNA", group.by = "patient", pt.size = 0.1) + 
  geom_hline(yintercept = c(corrected_count_threshold_low, corrected_count_threshold_high), 
             linetype = "dashed", color = "red") +
  ggtitle("Counts per Cell by Patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Figure 1B: Feature distribution
p2 <- VlnPlot(seurat, features = "nFeature_RNA", group.by = "patient", pt.size = 0.1) + 
  geom_hline(yintercept = corrected_feature_threshold_low, linetype = "dashed", color = "red") +
  ggtitle("Genes per Cell by Patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Figure 1C: Complexity distribution
p3 <- VlnPlot(seurat, features = "complexity", group.by = "patient", pt.size = 0.1) + 
  geom_hline(yintercept = corrected_complexity_threshold, linetype = "dashed", color = "red") +
  ggtitle("Complexity Score by Patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Figure 1D: Feature scatter with filtering boundaries
p4 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                     group.by = "patient") +
  geom_vline(xintercept = c(corrected_count_threshold_low, corrected_count_threshold_high), 
             linetype = "dashed", color = "red") +
  geom_hline(yintercept = corrected_feature_threshold_low, linetype = "dashed", color = "red") +
  ggtitle("Counts vs Genes with Filter Boundaries")

# Combine plots
figure1 <- (p1 | p2 | p3) / p4 + 
  plot_annotation(tag_levels = 'A', 
                  title = 'Figure 1: Quality Control Metrics with Filtering Thresholds')

ggsave("Figure1_QC_Metrics_Corrected.png", figure1, width = 16, height = 12)
print(figure1)



# Verify filtered data looks good
cat("\n=== FILTERED DATA SUMMARY ===\n")
cat("Cells after filtering:", ncol(seurat_filtered), "\n")
cat("Genes after filtering:", nrow(seurat_filtered), "\n")

cat("\nFiltered distributions:\n")
cat("nCount_RNA range:", range(seurat_filtered$nCount_RNA), "\n")
cat("nFeature_RNA range:", range(seurat_filtered$nFeature_RNA), "\n")
cat("Complexity range:", range(seurat_filtered$complexity), "\n")

# Check patient distribution after filtering
cat("\nPatient distribution after filtering:\n")
patient_dist <- table(seurat_filtered$patient)
print(patient_dist)

# Check what was removed
cells_removed <- setdiff(colnames(seurat), colnames(seurat_filtered))
cat("\nCells removed:", length(cells_removed), "\n")

# Check reasons for removal
removed_df <- data.frame(
  cell = cells_removed,
  low_counts = seurat$nCount_RNA[cells_removed] < corrected_count_threshold_low,
  high_counts = seurat$nCount_RNA[cells_removed] > corrected_count_threshold_high,
  low_features = seurat$nFeature_RNA[cells_removed] < corrected_feature_threshold_low,
  low_complexity = seurat$complexity[cells_removed] < corrected_complexity_threshold
)

cat("\nReasons for cell removal:\n")
print(colSums(removed_df[, -1]))

# Proceed with SCTransform normalization
seurat_filtered <- SCTransform(seurat_filtered, verbose = FALSE)

cat("Normalization completed successfully!\n")
cat("Variable features identified:", length(VariableFeatures(seurat_filtered)), "\n")

# Continue with PCA
seurat_filtered <- RunPCA(seurat_filtered, npcs = 50)
cat("PCA completed. Ready for clustering analysis.\n")


# Determine optimal dimensionality
pca_elbow <- ElbowPlot(seurat_filtered, ndims = 50) + 
  ggtitle("Figure 2A: PCA Elbow Plot - Dimensionality Assessment")
print(pca_elbow)

# Check PCA heatmaps for biological interpretation
pca_heatmap <- DimHeatmap(seurat_filtered, dims = 1:12, cells = 500, balanced = TRUE)
print(pca_heatmap)

# Proceed with clustering using 20-30 PCs based on elbow point
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:25)
seurat_filtered <- FindClusters(seurat_filtered, resolution = c(0.2, 0.5, 0.8, 1.2))

# Set the main clustering resolution (start with 0.5)
Idents(seurat_filtered) <- "SCT_snn_res.0.5"
cat("Clustering completed. Number of clusters:", length(unique(Idents(seurat_filtered))), "\n")

# Run UMAP
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:25)


# Figure 2A: UMAP colored by clusters
p2a <- DimPlot(seurat_filtered, label = TRUE, repel = TRUE) + 
  ggtitle("Cell Clusters") +
  theme(legend.position = "right")

# Figure 2B: UMAP colored by patient
p2b <- DimPlot(seurat_filtered, group.by = "patient") + 
  ggtitle("Colored by Patient") +
  theme(legend.position = "right")

# Figure 2C: UMAP colored by experiment
p2c <- DimPlot(seurat_filtered, group.by = "experiment") + 
  ggtitle("Colored by Experiment") +
  theme(legend.position = "right")

# Figure 2D: UMAP colored by depletion status
p2d <- DimPlot(seurat_filtered, group.by = "depletion_batch") + 
  ggtitle("Colored by Depletion Status") +
  theme(legend.position = "right")

# Combine Figure 2
figure2 <- (p2a | p2b) / (p2c | p2d) + 
  plot_annotation(tag_levels = 'A', 
                  title = 'Figure 2: Cell Clustering and Batch Effects')

ggsave("Figure2_Clustering_UMAP.png", figure2, width = 14, height = 12)
print(figure2)


# Find marker genes for each cluster
cluster_markers <- FindAllMarkers(seurat_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cat("Top markers per cluster:\n")
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
print(top_markers)

# Create dot plot of known cell type markers
epithelial_markers <- c("KRT19", "KRT8", "KRT18", "EPCAM", "CDH1")
immune_markers <- c("PTPRC", "CD3D", "CD79A", "CD68", "CD14") 
stromal_markers <- c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2")

# Figure 3A: Dot plot of cell type markers
p3a <- DotPlot(seurat_filtered, features = c(epithelial_markers, immune_markers, stromal_markers)) +
  RotatedAxis() + 
  ggtitle("Cell Type Marker Expression")

print(p3a)
