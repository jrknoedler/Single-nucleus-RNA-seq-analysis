#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_PrimedvUnprimed_Exploratory/POA_PrimedUnprimed_Exploratory_minimal_filtered"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_PrimedvUnprimed_Exploratory/POA_PrimedUnprimed_Exploratory_minimal.rds")
mySeurat
filteredSeurat <- subset(mySeurat, idents=c("9", "21"), invert=TRUE)
filteredSeurat <- RunPCA(filteredSeurat, verbose=TRUE)
filteredSeurat <- FindNeighbors(filteredSeurat, dims=1:30, verbose=TRUE)
filteredSeurat <- FindClusters(filteredSeurat, verbose=TRUE, resolution=0.5)
filteredSeurat <- RunUMAP(filteredSeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE, split.by="state")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE, group.by="state")
dev.off()
DefaultAssay(filteredSeurat) <- "RNA"
filteredSeurat.markers <- FindAllMarkers(filteredSeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(filteredSeurat.markers, file=paste0(output, "_allposmarkers.csv"))
DefaultAssay(filteredSeurat) <- "SCT"
pdf(paste0(output, "_TopMarkerheatmap.pdf"))
top10 <- filteredSeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(filteredSeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(filteredSeurat, file=(paste0(output, ".rds")))
