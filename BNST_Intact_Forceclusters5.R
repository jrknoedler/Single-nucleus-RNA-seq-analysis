#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_40pcs_res1_lowthresh_res1.2"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

Intact <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1__Intact.rds")
Intact <- SCTransform(Intact, verbose=TRUE)
Intact <- RunPCA(Intact)
warnings()
Intact <- FindNeighbors(Intact, dims=1:40)
Intact <- FindClusters(Intact, resolution=1.2)

Intact <- RunUMAP(Intact, dims=1:40)
pdf(paste0(output,"_Intact_UMAP.pdf"))
DimPlot(Intact, reduction="umap", label=TRUE)
dev.off()
Intact <- BuildClusterTree(Intact, dims=1:40)
pdf(paste0(output,"_Intact_clustertree.pdf"))
PlotClusterTree(object=Intact)
dev.off()
Intact.markers <- FindAllMarkers(Intact, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.2, test.use="MAST")
write.csv(Intact.markers, file=paste0(output,"_Intact_allposmarkers.csv"))
pdf(paste0(output,"_Intact_TopMarkerheatmap.pdf"))
top10 <- Intact.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Intact, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Intact, file=(paste0(output, "_Intact.rds")))
