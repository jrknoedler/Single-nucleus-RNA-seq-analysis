#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

Intact <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentRound1__Intact.rds")
Primed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentRound1__Primed.rds")
Unprimed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentRound1__Unprimed.rds")
Intact <- subset(Intact, idents=c("0","21"), invert=TRUE)
Primed <- subset(Primed, idents=c("0"), invert=TRUE)
Unprimed <- subset(Unprimed, idents=c("0"), invert=TRUE)

#Retransform, recluster and save sample 1
Primed <- SCTransform(Primed, verbose=TRUE)
Primed <- RunPCA(Primed)
warnings()
Primed <- FindNeighbors(Primed, dims=1:30)
Primed <- FindClusters(Primed, resolution=1)
Primed <- RunUMAP(Primed, dims=1:30)
pdf(paste0(output,"_Primed_UMAP.pdf"))
DimPlot(Primed, reduction="umap", label=TRUE)
dev.off()
Primed.markers <- FindAllMarkers(Primed, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Primed.markers, file=paste0(output,"_Primed_allposmarkers.csv"))
pdf(paste0(output,"_Primed_TopMarkerheatmap.pdf"))
top10 <- Primed.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Primed, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Primed, file=(paste0(output, "_Primed.rds")))

#Retransform, recluster and save sample 2
Unprimed <- SCTransform(Unprimed, verbose=TRUE)
Unprimed <- RunPCA(Unprimed)
warnings()
Unprimed <- FindNeighbors(Unprimed, dims=1:30)
Unprimed <- FindClusters(Unprimed, resolution=1)
Unprimed <- RunUMAP(Unprimed, dims=1:30)
pdf(paste0(output,"_Unprimed_UMAP.pdf"))
DimPlot(Unprimed, reduction="umap", label=TRUE)
dev.off()
Unprimed.markers <- FindAllMarkers(Unprimed, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Unprimed.markers, file=paste0(output,"_Unprimed_allposmarkers.csv"))
pdf(paste0(output,"_Unprimed_TopMarkerheatmap.pdf"))
top10 <- Unprimed.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Unprimed, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Unprimed, file=(paste0(output, "_Unprimed.rds")))

#Retransform, cluster and save sample3
Intact <- SCTransform(Intact, verbose=TRUE)
Intact <- RunPCA(Intact)
warnings()
Intact <- FindNeighbors(Intact, dims=1:30)
Intact <- FindClusters(Intact, resolution=1)
Intact <- RunUMAP(Intact, dims=1:30)
pdf(paste0(output,"_Intact_UMAP.pdf"))
DimPlot(Intact, reduction="umap", label=TRUE)
dev.off()
Intact.markers <- FindAllMarkers(Intact, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Intact.markers, file=paste0(output,"_Intact_allposmarkers.csv"))
pdf(paste0(output,"_Intact_TopMarkerheatmap.pdf"))
top10 <- Intact.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Intact, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Intact, file=(paste0(output, "_Intact.rds")))



