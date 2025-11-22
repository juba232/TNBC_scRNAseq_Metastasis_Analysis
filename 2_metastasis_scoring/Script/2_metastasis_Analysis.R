# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
library(reshape2)

# Load the processed Seurat object
seurat_filtered <- readRDS("seurat_processed_metastasis.rds")

# Improved gene set loading function
load_metastasis_genes <- function(file_path) {
  genes <- readLines(file_path)
  genes <- genes[genes != "Gene"]  # Remove header
  genes <- trimws(genes)  # Remove whitespace
  genes <- genes[genes != ""]  # Remove empty lines
  genes <- genes[!is.na(genes)]  # Remove NAs
  return(genes)
}


### YOU MUST HAVE THE THE METASTASIS FOLDER IN THE DIRECTORY ##


# Load gene sets with better error handling
artega_sig <- load_metastasis_genes("metastasis_genesets/artega_sig.txt")
mammaprint_sig <- load_metastasis_genes("metastasis_genesets/mammaprint_sig_new.txt") 
werb_sig <- load_metastasis_genes("metastasis_genesets/werb_49_metastasis_sig.txt")

cat("Gene sets loaded:\n")
cat("Artega genes:", length(artega_sig), "\n")
cat("MammaPrint genes:", length(mammaprint_sig), "\n")
cat("Werb genes:", length(werb_sig), "\n")

# Check which genes are actually present in our dataset
check_gene_presence <- function(gene_list, object_genes) {
  present_genes <- gene_list[gene_list %in% object_genes]
  missing_genes <- gene_list[!gene_list %in% object_genes]
  
  cat("Genes present:", length(present_genes), "/", length(gene_list), "\n")
  if(length(missing_genes) > 0) {
    cat("First 10 missing genes:", head(missing_genes, 10), "\n")
  }
  return(present_genes)
}

cat("\nChecking gene presence in dataset:\n")
cat("=== Artega Signature ===\n")
artega_present <- check_gene_presence(artega_sig, rownames(seurat_filtered))

cat("\n=== MammaPrint Signature ===\n")
mammaprint_present <- check_gene_presence(mammaprint_sig, rownames(seurat_filtered))

cat("\n=== Werb Signature ===\n")
werb_present <- check_gene_presence(werb_sig, rownames(seurat_filtered))


# Calculate scores only with genes that are present
seurat_filtered <- AddModuleScore(
  seurat_filtered,
  features = list(
    Artega_Metastasis = artega_present,
    MammaPrint_Metastasis = mammaprint_present,
    Werb_Metastasis = werb_present
  ),
  name = c("Artega", "MammaPrint", "Werb")
)

cat("Metastasis scores calculated successfully!\n")

# Check the score columns were added
cat("New metadata columns:\n")
print(tail(colnames(seurat_filtered@meta.data), 10))

# Find the exact column names for metastasis scores
metastasis_cols <- grep("Artega|MammaPrint|Werb", colnames(seurat_filtered@meta.data), value = TRUE)
cat("\nMetastasis score columns found:\n")
print(metastasis_cols)

# Check the actual column names that were created
score_cols <- tail(colnames(seurat_filtered@meta.data), 3)
cat("\nLast 3 columns (likely the scores):\n")
print(score_cols)

# Let's examine these columns to confirm they contain the scores
cat("\nFirst few values of potential score columns:\n")
for(col in score_cols) {
  cat(col, "range:", range(seurat_filtered@meta.data[[col]]), "\n")
}


# Figure 3: UMAP with metastasis scores
p_artega <- FeaturePlot(seurat_filtered, features = "Artega1") + 
  ggtitle("Artega Metastasis Score") +
  scale_color_viridis_c(option = "plasma") +
  theme(legend.position = "bottom")

p_mammaprint <- FeaturePlot(seurat_filtered, features = "MammaPrint2") + 
  ggtitle("MammaPrint Metastasis Score") +
  scale_color_viridis_c(option = "plasma") +
  theme(legend.position = "bottom")

p_werb <- FeaturePlot(seurat_filtered, features = "Werb3") + 
  ggtitle("Werb Metastasis Score") +
  scale_color_viridis_c(option = "plasma") +
  theme(legend.position = "bottom")

# Combine UMAP plots
figure3 <- p_artega | p_mammaprint | p_werb
figure3 <- figure3 + plot_annotation(tag_levels = 'A', 
                                     title = 'Figure 3: Metastasis Signature Scores')

#ggsave("Figure3_Metastasis_Scores_UMAP.png", figure3, width = 18, height = 6)
print(figure3)

# Figure 4: Violin plots by cluster
v1 <- VlnPlot(seurat_filtered, features = "Artega1", pt.size = 0.1) + 
  ggtitle("Artega Score") + NoLegend() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
v2 <- VlnPlot(seurat_filtered, features = "MammaPrint2", pt.size = 0.1) + 
  ggtitle("MammaPrint Score") + NoLegend() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
v3 <- VlnPlot(seurat_filtered, features = "Werb3", pt.size = 0.1) + 
  ggtitle("Werb Score") + NoLegend() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure4 <- v1 | v2 | v3
figure4 <- figure4 + plot_annotation(tag_levels = 'A', 
                                     title = 'Figure 4: Metastasis Scores by Cluster')

#ggsave("Figure4_Metastasis_Violin.png", figure4, width = 15, height = 5)
print(figure4)


# Calculate and visualize average scores per cluster
cluster_assignments <- Idents(seurat_filtered)
metastasis_scores <- seurat_filtered@meta.data[, c("Artega1", "MammaPrint2", "Werb3")]

avg_scores <- data.frame(
  cluster = levels(cluster_assignments),
  Artega = tapply(metastasis_scores$Artega1, cluster_assignments, mean),
  MammaPrint = tapply(metastasis_scores$MammaPrint2, cluster_assignments, mean),
  Werb = tapply(metastasis_scores$Werb3, cluster_assignments, mean)
)

# Add combined score
avg_scores$Combined <- rowMeans(avg_scores[, c("Artega", "MammaPrint", "Werb")])

cat("Average metastasis scores by cluster:\n")
print(avg_scores)

# Identify high-metastasis clusters (top 25% for each signature)
high_artega <- avg_scores$cluster[avg_scores$Artega > quantile(avg_scores$Artega, 0.75)]
high_mammaprint <- avg_scores$cluster[avg_scores$MammaPrint > quantile(avg_scores$MammaPrint, 0.75)]
high_werb <- avg_scores$cluster[avg_scores$Werb > quantile(avg_scores$Werb, 0.75)]
high_combined <- avg_scores$cluster[avg_scores$Combined > quantile(avg_scores$Combined, 0.75)]

cat("\nHigh metastasis clusters:\n")
cat("Artega high:", paste(high_artega, collapse = ", "), "\n")
cat("MammaPrint high:", paste(high_mammaprint, collapse = ", "), "\n")
cat("Werb high:", paste(high_werb, collapse = ", "), "\n")
cat("Combined high:", paste(high_combined, collapse = ", "), "\n")



# Figure 5: Heatmap of metastasis scores by cluster
avg_scores_melted <- melt(avg_scores, id.vars = "cluster", 
                          variable.name = "Signature", value.name = "Score")

figure5 <- ggplot(avg_scores_melted, aes(x = cluster, y = Signature, fill = Score)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Score") +
  theme_minimal() +
  labs(title = "Figure 5: Metastasis Potential by Cluster",
       x = "Cluster", y = "Metastasis Signature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_text(aes(label = round(Score, 2)), color = "white", size = 3)

#ggsave("Figure5_Metastasis_Heatmap.png", figure5, width = 10, height = 6)
print(figure5)


# Define high metastasis cells based on combined score
seurat_filtered$combined_metastasis <- rowMeans(metastasis_scores)
high_metastasis_threshold <- quantile(seurat_filtered$combined_metastasis, 0.75)
seurat_filtered$metastasis_status <- ifelse(seurat_filtered$combined_metastasis > high_metastasis_threshold, 
                                            "High", "Low")

# Figure 6: Metastatic subpopulations
p_status <- DimPlot(seurat_filtered, group.by = "metastasis_status", cols = c("blue", "red")) +
  ggtitle("Figure 6: Metastatic Subpopulations") +
  labs(subtitle = paste("High metastasis:", sum(seurat_filtered$metastasis_status == "High"), "cells"))

# Add proportion by cluster
cluster_metastasis <- table(Idents(seurat_filtered), seurat_filtered$metastasis_status)
cluster_metastasis_prop <- prop.table(cluster_metastasis, margin = 1)
high_metastasis_prop <- cluster_metastasis_prop[, "High"]

prop_df <- data.frame(
  cluster = names(high_metastasis_prop),
  proportion = high_metastasis_prop
)

p_prop <- ggplot(prop_df, aes(x = cluster, y = proportion, fill = proportion)) +
  geom_col() +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  labs(title = "Proportion of High-Metastasis Cells",
       x = "Cluster", y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure6 <- p_status | p_prop
figure6 <- figure6 + plot_annotation(title = 'Figure 6: Metastatic Subpopulation Analysis')
#ggsave("Figure6_Metastasis_Subpopulations.png", figure6, width = 14, height = 6)
print(figure6)


# Figure 7: Patient-specific metastasis patterns
# Average metastasis scores by patient
patient_scores <- seurat_filtered@meta.data %>%
  group_by(patient) %>%
  summarise(
    Artega = mean(Artega1),
    MammaPrint = mean(MammaPrint2),
    Werb = mean(Werb3),
    Combined = mean(combined_metastasis),
    n_cells = n()
  )

cat("\nAverage metastasis scores by patient:\n")
print(patient_scores)

# Visualize patient patterns
p_patient <- ggplot(patient_scores, aes(x = patient, y = Combined, fill = patient)) +
  geom_col() +
  labs(title = "Figure 7: Average Metastasis Potential by Patient",
       x = "Patient", y = "Combined Metastasis Score") +
  theme_minimal() +
  theme(legend.position = "none")

# Patient distribution in high vs low metastasis
patient_metastasis <- table(seurat_filtered$patient, seurat_filtered$metastasis_status)
patient_metastasis_prop <- prop.table(patient_metastasis, margin = 1)

prop_melted <- melt(patient_metastasis_prop, varnames = c("patient", "status"), value.name = "proportion")

p_patient_prop <- ggplot(prop_melted, aes(x = patient, y = proportion, fill = status)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("Low" = "blue", "High" = "red")) +
  labs(title = "Metastasis Status Distribution by Patient",
       x = "Patient", y = "Proportion") +
  theme_minimal()

figure7 <- p_patient | p_patient_prop
figure7 <- figure7 + plot_annotation(title = 'Figure 7: Patient-Specific Metastasis Patterns')
#ggsave("Figure7_Patient_Metastasis.png", figure7, width = 14, height = 6)
print(figure7)



# Generate comprehensive summary
cat("\n=== TRIPLE-NEGATIVE BREAST CANCER METASTASIS ANALYSIS SUMMARY ===\n")
cat("Dataset: ", ncol(seurat_filtered), " cells from ", length(unique(seurat_filtered$patient)), " patients\n")
cat("Clusters identified: ", length(unique(Idents(seurat_filtered))), "\n")
cat("High-metastasis cells: ", sum(seurat_filtered$metastasis_status == "High"), 
    " (", round(mean(seurat_filtered$metastasis_status == "High") * 100, 1), "%)\n", sep = "")

cat("\nKey Findings:\n")
cat("• High-metastasis clusters: ", paste(high_combined, collapse = ", "), "\n")
cat("• Patient with highest metastasis potential: ", 
    patient_scores$patient[which.max(patient_scores$Combined)], "\n")
cat("• Most consistent high-metastasis signature: ", 
    names(which.max(c(Artega = length(high_artega), MammaPrint = length(high_mammaprint), Werb = length(high_werb)))), "\n")

cat("\nFigures Generated:\n")
cat("• Figure 3: Metastasis signature scores on UMAP\n")
cat("• Figure 4: Violin plots of scores by cluster\n")
cat("• Figure 5: Heatmap of cluster-level metastasis potential\n")
cat("• Figure 6: Metastatic subpopulation identification\n")
cat("• Figure 7: Patient-specific metastasis patterns\n")

# Save final object with all analysis
saveRDS(seurat_filtered, "TNBC_metastasis_analysis_complete.rds")
cat("\nComplete analysis saved to: TNBC_metastasis_analysis_complete.rds\n")
