#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/UnprimedBNST_Preprocessing/UnprimedBNST_hemo"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/UnprimedBNST_Preprocessing/UnprimedBNST_preprocessing.rds")
mySeurat
filteredSeurat <- subset(mySeurat, idents=c("2", "22","23"), invert=TRUE)
filteredSeurat
filteredSeurat <- RunPCA(filteredSeurat, verbose=TRUE)
filteredSeurat <- FindNeighbors(filteredSeurat, dims=1:30, verbose=TRUE)
filteredSeurat <- FindClusters(filteredSeurat, verbose=TRUE, resolution=0.5)
filteredSeurat <- RunUMAP(filteredSeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(filteredSeurat, reduction="umap", do.label=TRUE)
dev.off()
DefaultAssay(filteredSeurat) <- "RNA"
NormalizeData(filteredSeurat)
pdf(paste0(output,"_Hbb-bsfeatureplot.pdf"))
FeaturePlot(filteredSeurat, features=c("Hbb-bs"))
dev.off()
pdf(paste0(output, "_Hbb-bsvln.pdf"))
VlnPlot(filteredSeurat, features =c("Hbb-bs"), pt.size=0)
dev.off()
pdf(paste0(output,"_Hba-a1featureplot.pdf"))
FeaturePlot(filteredSeurat, features=c("Hba-a1"))
dev.off()
pdf(paste0(output, "_Hba-a1svln.pdf"))
VlnPlot(filteredSeurat, features =c("Hba-a1"), pt.size=0)
dev.off()
pdf(paste0(output,"_Hba-a2featureplot.pdf"))
FeaturePlot(filteredSeurat, features=c("Hba-a2"))
dev.off()
pdf(paste0(output, "_Hba-a2svln.pdf"))
VlnPlot(filteredSeurat, features =c("Hba-a2"), pt.size=0)
dev.off()
