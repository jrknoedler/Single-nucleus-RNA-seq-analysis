#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_Primed_40pcs_res1.5_lowthresh"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

Primed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1__Primed.rds")
Primed <- SCTransform(Primed, verbose=TRUE)
Primed <- RunPCA(Primed)
warnings()
Primed <- FindNeighbors(Primed, dims=1:40)
Primed <- FindClusters(Primed, resolution=1.5)

Primed <- RunUMAP(Primed, dims=1:40)
pdf(paste0(output,"_Primed_UMAP.pdf"))
DimPlot(Primed, reduction="umap", label=TRUE)
dev.off()
Primed <- BuildClusterTree(Primed, dims=1:40)
pdf(paste0(output,"_Primed_clustertree.pdf"))
PlotClusterTree(object=Primed)
dev.off()
Primed.markers <- FindAllMarkers(Primed, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.20, test.use="MAST")
write.csv(Primed.markers, file=paste0(output,"_Primed_allposmarkers.csv"))
pdf(paste0(output,"_Primed_TopMarkerheatmap.pdf"))
top10 <- Primed.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Primed, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Primed, file=(paste0(output, "_Primed.rds")))
