#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/MeA_DEGexclude_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")

DefaultAssay(mySeurat) <- "integrated"
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)

genelist2 <- read.table("topGO/Total_ByRegion/MeA_Genesonly.txt", header=FALSE)
genelist2 <- unlist(genelist2)
genelist2 <- as.matrix(genelist2)

hvg <- mySeurat@assays$integrated@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.int <- setdiff(hvg, genelist)
hvg.final <- setdiff(hvg.int, genelist2)
mySeurat <- RunPCA(mySeurat, features=hvg.final)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
saveRDS(mySeurat, file=paste0(output, ".rds"))