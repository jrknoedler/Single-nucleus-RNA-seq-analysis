#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/Sakura_IntegrationTests/Naive_Integration_Fixed/BNST/IntegratedBNST_Naive_Filtered1"

library(Seurat)
library(rgeos)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)


mySeurat <- readRDS("Seurat/Sakura_IntegrationTests/Naive_Integration_Fixed/BNST/BNST_Naive_Integration_PC30_Res12_Fixed.rds")
mySeurat <- subset(mySeurat, idents=c("26","37","42","43"), invert=TRUE)
DefaultAssay(mySeurat) <- "SCT"
mySeurat <- RunPCA(mySeurat)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.2)
saveRDS(mySeurat, file=paste0(output, ".rds"))
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=40)
VlnPlot(mySeurat, features=c("Slc23a1","Slc17a6","Esr1"), ncol=1, pt.size=0)
dev.off()
SeuratScaled <- ScaleData(mySeurat)
#pdf(paste0(output,"_Markerdot.pdf"), width=40, height=40)
#DotPlot(SeuratScaled, features=c("Gad1","Slc17a6","Esr1","Cckar","Mbp","Plp1","Cxcrl1","Inpp5d"))
#dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output,"RNAfeat.pdf"), height=7, width=42)
VlnPlot(mySeurat, features=("nCount_RNA"), pt.size=0)
dev.off()
pdf(paste0(output,"genenum.pdf"), height=7, width=42)
VlnPlot(mySeurat, features=c("nFeature_RNA"), pt.size=0)
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.40, logfc.threshold = 0.40, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
