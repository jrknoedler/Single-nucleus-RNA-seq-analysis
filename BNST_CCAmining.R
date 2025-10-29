#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNSTMacosko/CCA_1stpassanalysisres0.8"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

BNSTCCA <- readRDS("Seurat/BNSTMacosko/Macosko_JRK_Integrated_firstpass.rds")
BNSTCCA <- RunUMAP(BNSTCCA, reduction="pca", dims=1:30)
BNSTCCA <- FindNeighbors(BNSTCCA, reduction ="pca", dims=1:30)
BNSTCCA <- FindClusters(BNSTCCA, resolution=0.8)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(BNSTCCA, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_modesplit.pdf"))
DimPlot(BNSTCCA, reduction="umap", label=TRUE, split.by="type")
dev.off()
pdf(paste0(output,"_UMAP_modelabel.pdf"))
DimPlot(BNSTCCA, reduction="umap", label=TRUE, group.by="type")
dev.off()
DefaultAssay(BNSTCCA) <- "RNA"
#BNSTCCA2  <- ScaleData(BNSTCCA)
#saveRDS(BNSTCCA2, file=paste0(output,"zscored.rds"))
Markers <- read.table("DotPlotMarkerLists/BNST_ReorderedFinalalternate_CPM.txt", header=FALSE)
Markers <- unlist(Markers)
pdf(file=paste0(output,"Vlnplottest_res0.8.pdf"), width=30, height=200)
VlnPlot(BNSTCCA, features=Markers, split.by="type", ncol=1, pt.size=0)
dev.off()
