#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_PrimedvUnprimed_doubletsremoved/PrimedvUnprimed_doubletsremoved_filtered1"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_PrimedvUnprimed_doubletsremoved/POA_PrimedvUnprimed_Exploratory_doubletsremovedexploratory.rds")
mySeurat
filteredSeurat <- subset(mySeurat, idents=c("40","39","37"), invert=TRUE)
DefaultAssay(filteredSeurat) <- "integrated"
filteredSeurat <- RunPCA(filteredSeurat, verbose=TRUE)
filteredSeurat <- FindNeighbors(filteredSeurat, dims=1:30, verbose=TRUE)
filteredSeurat <- FindClusters(filteredSeurat, verbose=TRUE, resolution=1)
filteredSeurat <- RunUMAP(filteredSeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_statesplit.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE, split.by="state")
dev.off()
pdf(paste0(output,"_UMAP_statelabel.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE, group.by="state")
dev.off()
filteredSeurat.markers <- FindAllMarkers(filteredSeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(filteredSeurat.markers, file=paste0(output, "_allposmarkers.csv"))
saveRDS(filteredSeurat, file=(paste0(output, ".rds")))
