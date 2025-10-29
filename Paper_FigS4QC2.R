#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/PaperQC/PaperQC"
regoutput <- "RegulonsCompiled/POA/Sexmarkersperclust/"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)




mySeurat <- readRDS("Seurat/PaperQC/PaperQCprePCA.rds")
head(mySeurat[[]])
mySeurat[['seurat_clusters']] <- NULL
mySeurat[['SCT_snn_res.1']] <- NULL
mySeurat[['SCT_snn_res.1.2']] <- NULL
mySeurat[['SCT_snn_res.1.5']] <- NULL
mySeurat[['SCT_snn_res.0.4']] <- NULL
mySeurat[['SCT_snn_res.0.8']] <- NULL
mySeurat[['integrated_snn_res.1.2']] <- NULL
head(mySeurat[[]])
Idents(mySeurat) <- "sample"
my_levels <- c("MaleBNST","PrimedBNST","UnprimedBNST","MaleMeA","PrimedMeA","UnprimedMeA","IntactPOA","PrimedPOA","UnprimedPOA","MaleVMH","PrimedVMH","UnprimedVMH")
mySeurat@active.ident <- factor(x=mySeurat@active.ident, levels=my_levels)
pdf(file=paste0(output,"QCUMIplot.pdf"),width=15)
VlnPlot(mySeurat, features=c("nCount_RNA"),pt.size=0)
dev.off()
pdf(file=paste0(output,"QCgeneplot.pdf"),width=15)
VlnPlot(mySeurat, features=c("nFeature_RNA"), pt.size=0)
dev.off()
mySeurat <- RunPCA(mySeurat, dims=1:100)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dims=1:50)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:50)
mySeurat <- FindClusters(mySeurat, resolution=2)
pdf(paste0(output,"_UMAP.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()

pdf(paste0(output,"_UMAP_sexlabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", shuffle=TRUE, label=TRUE, group.by="sex", label=FALSE)
dev.off()

pdf(paste0(output,"_UMAP_hormonelabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", shuffle=TRUE, label=TRUE, group.by="Hormone", label=FALSE)
dev.off()


pdf(paste0(output,"_UMAP_regionlabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", shuffle=TRUE, label=TRUE, group.by="region", label=FALSE)
dev.off()

pdf(paste0(output,"_UMAP_batchlabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", shuffle=TRUE, label=TRUE, group.by="batch", label=FALSE)
dev.off()


pdf(paste0(output,"_UMAP_idabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sample", label=FALSE)
dev.off()
saveRDS(mySeurat, file=paste0(output, ".rds"))