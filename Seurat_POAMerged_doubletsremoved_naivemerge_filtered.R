#!/usr/bin/env Rscript

output <- "Seurat/POA_MalevUnprimedMerged_Exploratory/POA_MalevUnprimed_Exploratory_doubleltsremoved_naivemergefiltered_sexremove"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/POA_MalevUnprimed_Exploratory_doubleltsremoved_naivemerge.rds")
mySeurat <- subset(mySeurat, idents=c("0","24","32"), invert=TRUE)
mySeurat <- FindVariableFeatures(mySeurat, selection.method="vst", nfeatures=2000)
mySeurat <- ScaleData(mySeurat)
genelist <- read.table("Genelists/Sexchrom.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
mySeurat <- FindVariableFeatures(mySeurat, selection.method="vst", nfeatures=2000)
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$RNA@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
mySeurat <- RunPCA(mySeurat, features=hvg.final)
pdf(paste0(output,"_PCAplot.pdf"))
DimPlot(mySeurat, reduction="pca", group.by="sex")
dev.off()
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_origlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sex")
dev.off()
mySeuratSeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeuratSeurat.markers, file=paste0(output, "_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
