# 1_TNBC scRNA-seq Preprocessing
This folder contains QC, filtering, normalization, PCA, UMAP, and clustering steps.

## Workflow

- [x] Load raw count matrix
- [x] Quality control (counts, genes, % mito)
- [x] Filter cells & genes
- [x] Normalize using SCTransform
- [x] Run PCA (25 PCs)
- [x] Run UMAP
- [x] Identify clusters (resolution = 0.5)

