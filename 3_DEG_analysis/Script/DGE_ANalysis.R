# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)
library(reshape2)
library(EnhancedVolcano)



seurat_filtered <- load("TNBC_metastasis_analysis_complete.rds")

# Perform differential expression between high vs low metastasis cells
cat("Performing differential expression analysis...\n")
Idents(seurat_filtered) <- "metastasis_status"
de_metastasis <- FindMarkers(seurat_filtered, 
                             ident.1 = "High", 
                             ident.2 = "Low",
                             min.pct = 0.1,
                             logfc.threshold = 0.25,
                             test.use = "wilcox")

# Filter and sort significant results
de_metastasis_sig <- de_metastasis[de_metastasis$p_val_adj < 0.05, ]
de_metastasis_sig <- de_metastasis_sig[order(-de_metastasis_sig$avg_log2FC), ]

cat("Differential expression results:\n")
cat("Total significant genes (adj. p < 0.05):", nrow(de_metastasis_sig), "\n")
cat("Upregulated in high-metastasis:", sum(de_metastasis_sig$avg_log2FC > 0), "\n")
cat("Downregulated in high-metastasis:", sum(de_metastasis_sig$avg_log2FC < 0), "\n")

# Top upregulated genes in high-metastasis cells
cat("\nTop 20 upregulated genes in high-metastasis cells:\n")
top_upregulated <- head(de_metastasis_sig[de_metastasis_sig$avg_log2FC > 0, ], 20)
print(top_upregulated[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")])

# Top downregulated genes
cat("\nTop 20 downregulated genes in high-metastasis cells:\n")
top_downregulated <- head(de_metastasis_sig[de_metastasis_sig$avg_log2FC < 0, ], 20)
print(top_downregulated[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")])



# Figure 8: Volcano plot of differential expression
figure8 <- EnhancedVolcano(de_metastasis,
                           lab = rownames(de_metastasis),
                           x = 'avg_log2FC',
                           y = 'p_val_adj',
                           title = 'Figure 8: High vs Low Metastasis Potential - Differential Expression',
                           pCutoff = 0.05,
                           FCcutoff = 0.5,
                           pointSize = 2.0,
                           labSize = 4.0,
                           colAlpha = 0.7,
                           legendPosition = 'right',
                           drawConnectors = TRUE,
                           widthConnectors = 0.5)

#ggsave("Figure8_DE_Volcano.png", figure8, width = 12, height = 8)
print(figure8)

# Figure 9: Heatmap of top differential genes
top_genes <- rownames(head(de_metastasis_sig, 30))
Idents(seurat_filtered) <- "metastasis_status"

figure9 <- DoHeatmap(seurat_filtered, 
                     features = top_genes,
                     group.by = "metastasis_status",
                     group.colors = c("Low" = "blue", "High" = "red")) +
  ggtitle("Figure 9: Top Differential Genes - High vs Low Metastasis") +
  theme(axis.text.y = element_text(size = 8))

#ggsave("Figure9_DE_Heatmap.png", figure9, width = 10, height = 12)
print(figure9)
