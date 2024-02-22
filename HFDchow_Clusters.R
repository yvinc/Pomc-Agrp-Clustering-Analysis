library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library("dbscan")

# In case there's any randomization during the analysis, set seed:
set.seed(3598)

# Set working directory:
setwd("Desktop/REX/data")
ARC_ME <- "GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt"
counts_ARC_ME <- data.table::fread(ARC_ME)
# ncol(counts_ARC_ME) #20921 barcodes, processed

# Gene names are in the first column, so we need to move them to rownames:
counts_ARC_ME <- counts_ARC_ME %>% 
  column_to_rownames("V1")

# Uncomment below to view row/column names (gene names):
# rownames(counts_ARC_ME)
# colnames(counts_ARC_ME)

# Subset the table into hfd and chow condition dataframes:
counts_hfd <- as.data.frame(c(select(counts_ARC_ME, contains("hfd"))))
counts_chow <- as.data.frame(c(select(counts_ARC_ME, contains("chow"))))
combined_HFD_chow <- as.data.frame(c(counts_hfd, counts_chow))
rownames(combined_HFD_chow) <- rownames(counts_ARC_ME)
rm(counts_ARC_ME)
# Preview dataframe:
View(combined_HFD_chow)

# Create Seurat Object from dataframe `combined_HFD_chow`
HFD_chow <- CreateSeuratObject(counts = combined_HFD_chow, project = "HFDchow",
                               min.cells = 3, min.features = 200)

# No splitting in this file

## Seurat PC Analysis / UMAP Workflow:
## Replace value produced from Seurat's NormalizeData() with the original count values
# 1. Create separate count matrix object for both groups from counts
counts_all <- HFD_chow[["RNA"]]$counts
# 2. Run NormalizeData() on HFD_chow.
HFD_chow_wkfl <- NormalizeData(HFD_chow)
# 3. Examine values produced by NormalizeData() in the `data.Chow` and `data.HFD` slots/layers.
HFD_chow_wkfl[["RNA"]]$data 
# 4. Replace normalized data matrix with the previous count matrix.
HFD_chow_wkfl[["RNA"]]$data <- counts_all

# 5. FindVariableFeatures() and ScaleData
HFD_chow_wkfl <- FindVariableFeatures(HFD_chow_wkfl)
HFD_chow_wkfl <- ScaleData(HFD_chow_wkfl, do.scale = FALSE, do.center = TRUE) # Set do.center = TRUE, but set do.scale = FALSE
# Scaled data matrix is stored in HFD_chow_wkfl[["RNA"]]$scale.data
# Uncomment below to examine scaled values:
# View(HFD_chow_wkfl[["RNA"]]$scale.data)

# 6. Run PCA and generate graphs
HFD_chow_wkfl <- RunPCA(HFD_chow_wkfl)
VizDimLoadings(HFD_chow_wkfl, dims = 1:2, reduction = "pca")
DimPlot(HFD_chow_wkfl, reduction = "pca") + NoLegend()
ElbowPlot(HFD_chow_wkfl)

# 7. Find Nearest Neighbors and Clustering
HFD_chow_wkfl <- FindNeighbors(HFD_chow_wkfl, dims = 1:25, reduction = "pca")
HFD_chow_wkfl <- FindClusters(HFD_chow_wkfl, resolution = 0.8, cluster.name = "HFD_chow_clusters")

# 8. Run UMAP and generate plot
HFD_chow_UMAP <- RunUMAP(HFD_chow_wkfl, dims = 1:25)
DimPlot(HFD_chow_UMAP, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
DimPlot(HFD_chow_UMAP, reduction = "umap", label = TRUE)

