library(Seurat)
library(tidyverse)
library(data.table)
library(patchwork)
library(devtools)
library(ggpubr)

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
HFD_chow_wkfl <- ScaleData(HFD_chow_wkfl, do.scale = F, do.center = T) # Set do.center = TRUE, but set do.scale = FALSE
# Scaled data matrix is stored in HFD_chow_wkfl[["RNA"]]$scale.data
# Uncomment below to examine scaled values:
# View(HFD_chow_wkfl[["RNA"]]$scale.data)

# 6. Run PCA and generate graphs
HFD_chow_wkfl <- RunPCA(HFD_chow_wkfl)
VizDimLoadings(HFD_chow_wkfl, dims = 1:2, reduction = "pca")
PCA_plot <- DimPlot(HFD_chow_wkfl, reduction = "pca") + ggtitle("Principle Component Analysis (PCA)") + 
              theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
PCA_plot
ggsave("PCA_plot.png", plot = PCA_plot, dpi = 300, height = 13, width = 16, unit = 'cm')

ElbowPlot(HFD_chow_wkfl)

# 7. Find Nearest Neighbors and Clustering
HFD_chow_wkfl <- FindNeighbors(HFD_chow_wkfl, dims = 1:25, reduction = "pca")
HFD_chow_wkfl <- FindClusters(HFD_chow_wkfl, resolution = 0.8, cluster.name = "HFD_chow_clusters")

# 8. Run UMAP and generate plot
HFD_chow_UMAP <- RunUMAP(HFD_chow_wkfl, dims = 1:25)
DimPlot(HFD_chow_UMAP, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"))
UMAP_no_lable <- DimPlot(HFD_chow_UMAP, reduction = "umap", label = TRUE) + 
                  ggtitle("UMAP Dimention Reduction Clustering (PC = 25)") + 
                  theme(plot.title = element_text(hjust = 0.5))
UMAP_no_lable
ggsave("UMAP_no_lable.png", plot = UMAP_no_lable, dpi = 300, height = 13, width = 16, unit = 'cm')


# 9. Find specific expression of Pomc and Agrp in among the clusters: Identify which cluster has the highest expression of POMC/Agrp
heatmap_pomc1 <- FeaturePlot(HFD_chow_UMAP, features = "Pomc")
heatmap_agrp1 <- FeaturePlot(HFD_chow_UMAP, features = "Agrp")

heatmap_UMAP <- ggarrange(heatmap_pomc1, heatmap_agrp1,
          ncol = 1, nrow = 2,
          common.legend = TRUE, legend = "right")

ggsave("heatmap_UMAP.png", plot = heatmap_UMAP, dpi = 300, height = 20, width = 13, unit = 'cm')
ggsave("heatmap_pomc1.png", plot = heatmap_pomc1, dpi = 300, height = 11, width = 15, unit = 'cm')
ggsave("heatmap_agrp1.png", plot = heatmap_agrp1, dpi = 300, height = 11, width = 15, unit = 'cm')

# 10. Assign/Rename cell type identity (i.e. Agrp and Pomc) to cluster 4 and cluster 10.
# Rename cluster Idents with Agrp and Pomc and Display previous UMAP with new Agrp and Pomc cluster identity.
new.cluster.ids <- c("0", "1", "2", "3", "Agrp",
                     "5", "6", "7", "8", "9",
                     "Pomc", "11", "12", "13", "14",
                     "15", "16", "17", "18", "19", "20")
names(new.cluster.ids) <- levels(HFD_chow_UMAP)
HFD_chow_UMAP <- RenameIdents(HFD_chow_UMAP, new.cluster.ids) 

UMAP_labled <- DimPlot(HFD_chow_UMAP, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Dimention Reduction Clustering (PC = 25)") + theme(plot.title = element_text(hjust = 0.5))
UMAP_labled
ggsave("UMAP_labled.png", plot = UMAP_labled, dpi = 300, height = 13, width = 16, unit = 'cm')

# 11. Subset Seurat object into Agrp and Pomc cell type clusters only. 
pomc_agrp_subset <- subset(x = HFD_chow_UMAP, idents = c("Agrp", "Pomc"))
  # Compare the strucutre of previous object vs. the subset: Subset data has only 754 rows and 5 cols!
View(HFD_chow_UMAP)
View(pomc_agrp_subset)
table(Idents(pomc_agrp_subset))
levels(pomc_agrp_subset)
head(Idents(pomc_agrp_subset))

# 12. Repeat Find Neighbors and Clusters within the Agrp and Pomc clusters
pomc_agrp_wkfl <- FindNeighbors(pomc_agrp_subset, dims = 1:25, reduction = "pca")
pomc_agrp_wkfl <- FindClusters(pomc_agrp_wkfl, resolution = 0.8, cluster.name = "unintegrated_clusters")

# 13. Run UMAP and heatmap
pomc_agrp_UMAP <- RunUMAP(pomc_agrp_wkfl, dims = 1:25)
subset_UMAP <- DimPlot(pomc_agrp_UMAP, reduction = "umap", label = TRUE)
subset_UMAP <- subset_UMAP + ggtitle("Subclustered UMAP for Pomc and Agrp") + 
  theme(plot.title = element_text(hjust = 0.5))
subset_UMAP
ggsave("subset_UMAP.png", plot = subset_UMAP, dpi = 300, height = 13, width = 14, unit = 'cm')

DimPlot(pomc_agrp_UMAP, reduction = "umap", split.by = "orig.ident", label = TRUE)
DimPlot(pomc_agrp_UMAP, reduction = "umap", group.by = c("orig.ident", "seurat_clusters"), label = TRUE)
heat_subset <- FeaturePlot(pomc_agrp_UMAP, features = c("Pomc", "Agrp"))
heat_subset
ggsave("heat_subset.png", plot = heat_subset, dpi = 300, height = 10, width = 16, unit = 'cm')
FeaturePlot(pomc_agrp_UMAP, features = "Lepr")
FeaturePlot(pomc_agrp_UMAP, features = c("Lepr", "Insr"), split.by = "orig.ident")
FeaturePlot(pomc_agrp_UMAP, features = c("Gad1", "Slc17a6"), split.by = "orig.ident")


# 14. Boxplot of gene expression of Lepr, Insr, Gad1, and Slc17a6 in the Agrp and Pomc cluster from iv) chow animals and v) hfd animals.
# Make sure key genes are still in the subset Pomc and Agrp Seurat object (pomc_agrp_subset):
   "Lepr" %in% Features(pomc_agrp_subset)
   "Insr" %in% Features(pomc_agrp_subset)
   "Gad1" %in% Features(pomc_agrp_subset)
"Slc17a6" %in% Features(pomc_agrp_subset)

# Truncate expression level to > 0.5
pomc_agrp_truncate <- subset(pomc_agrp_subset, subset = Lepr > 0.5 & Insr > 0.5 & Gad1 > 0.5 & Slc17a6 > 0.5)

# Produce violin plot for key genes in Agrp
Vln_lepr_agrp <- VlnPlot(pomc_agrp_truncate, features = "Lepr", idents = "Agrp", group.by = "orig.ident",
                         pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
Vln_insr_agrp <- VlnPlot(pomc_agrp_truncate, features = "Insr", idents = "Agrp", group.by = "orig.ident",
                         pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
Vln_gad1_agrp <- VlnPlot(pomc_agrp_truncate, features = "Gad1", idents = "Agrp", group.by = "orig.ident",
                         pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
Vln_Slc17a6_agrp <- VlnPlot(pomc_agrp_truncate, features = "Slc17a6", idents = "Agrp", group.by = "orig.ident",
                         pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
vln_agrp_4 <- ggarrange(Vln_lepr_agrp, Vln_insr_agrp, Vln_gad1_agrp, Vln_Slc17a6_agrp,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2,
                      common.legend = TRUE, legend = "right")
vln_agrp_4 <- annotate_figure(vln_agrp_4, top = text_grob("Key Gene Expression Level in Agrp Neuron", size = 20, face = "bold"))
ggsave("vln_agrp_4.png", plot = vln_agrp_4, dpi = 300, height = 20, width = 20, unit = 'cm')

# Produce violin plot for key genes in Pomc
Vln_lepr_pomc <- VlnPlot(pomc_agrp_truncate, features = "Lepr", idents = "Pomc", group.by = "orig.ident",
                         pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
Vln_insr_pomc <- VlnPlot(pomc_agrp_truncate, features = "Insr", idents = "Pomc", group.by = "orig.ident",
                         pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
Vln_gad1_pomc <- VlnPlot(pomc_agrp_truncate, features = "Gad1", idents = "Pomc", group.by = "orig.ident",
                         pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
Vln_Slc17a6_pomc <- VlnPlot(pomc_agrp_truncate, features = "Slc17a6", idents = "Pomc", group.by = "orig.ident",
                            pt.size = 0.1, cols=c("#9AD9FF", "#EAC0C0"), add.noise = FALSE) 
vln_pomc_4 <- ggarrange(Vln_lepr_pomc, Vln_insr_pomc, Vln_gad1_pomc, Vln_Slc17a6_pomc,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2,
                      common.legend = TRUE, legend = "right")
vln_pomc_4 <- annotate_figure(vln_pomc_4, top = text_grob("Key Gene Expression Level in Pomc Neuron", size = 20, face = "bold"))
ggsave("vln_pomc_4.png", plot = vln_pomc_4, dpi = 300, height = 20, width = 20, unit = 'cm')


# VlnPlot(pomc_agrp_truncate, features = c("Lepr", "Insr", "Gad1", "Slc17a6"), idents = "Agrp", 
#           group.by = "orig.ident", pt.size = 0.1, cols=c("#EAC0C0", "#9AD9FF"), add.noise = FALSE) 
# VlnPlot(pomc_agrp_truncate, features = c("Lepr", "Insr", "Gad1", "Slc17a6"), idents = "Pomc", 
#         group.by = "orig.ident", pt.size = 0.1, cols=c("#EAC0C0", "#9AD9FF"), add.noise = FALSE) 
# 

# pomc_agrp_subset.de.markers <- FindMarkers(pomc_agrp_subset, ident.1 = "Agrp", ident.2 = "Pomc")
# head(pomc_agrp_subset.de.markers, 14)
# pomc_agrp_subset.de.markers[c("Lepr", "Insr", "Gad1", "Slc17a6"),]



# FYI
pomc_agrp_subset$celltype.condition <- paste(Idents(pomc_agrp_subset), pomc_agrp_subset$orig.ident, sep = "_")
pomc_agrp_subset$celltype <- Idents(pomc_agrp_subset)
Idents(pomc_agrp_subset) <- "celltype.condition"
condition.diffgenes_agpr <- FindMarkers(pomc_agrp_subset, ident.1 = "Agrp_HFD", ident.2 = "Agrp_Chow", 
                                        logfc.threshold = 0.25, verbose = FALSE)
condition.diffgenes_pomc <- FindMarkers(pomc_agrp_subset, ident.1 = "Pomc_HFD", ident.2 = "Pomc_Chow", verbose = FALSE)

view(condition.diffgenes_agpr)
view(condition.diffgenes_pomc)

condition.diffgenes_agpr %>%
  filter(row.names(condition.diffgenes_agpr) %in% c('Lepr', 'Insr', 'Gad1', 'Slc17a6'))

condition.diffgenes_pomc %>%
  filter(row.names(condition.diffgenes_pomc) %in% c('Lepr', 'Insr', 'Gad1', 'Slc17a6'))


# write.csv(condition.diffgenes_agpr, "Agrp_Differential.csv", row.names= TRUE)
# write.csv(condition.diffgenes_pomc, "Agrp_Differential.csv", row.names= TRUE)

#---


