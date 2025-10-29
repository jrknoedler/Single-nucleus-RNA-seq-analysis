#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_Anderson/AndersonFemale_GliaFilt2"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)


mySeurat <- readRDS("Seurat/VMH_Anderson/AndersonFemale_GliaFilt1.rds")
mySeurat <- subset(mySeurat, idents=c("0","51","39","40","33","34","14","24","43","1","2","4","6","12","16","21","23","25","27","28","29","30","31","35","47","48","49","53","54","56"), invert=TRUE)
mySeurat[["percent.mt"]] <- PercentageFeatureSet(mySeurat, pattern = "mt-^")

mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.mt"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
mySeurat <- RunPCA(mySeurat, features=hvg.final)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:50)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:50)
mySeurat <- FindClusters(mySeurat, resolution=1.2)
saveRDS(mySeurat, file=paste0(output, ".rds"))
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=40)
VlnPlot(mySeurat, features=c("Slc23a1","Slc17a6","Esr1","Cckar"), ncol=1, pt.size=0)
dev.off()
SeuratScaled <- ScaleData(mySeurat)
pdf(paste0(output,"_Markerdot.pdf"), width=40, height=40)
DotPlot(SeuratScaled, features=c("Gad1","Slc17a6","Esr1","Cckar","Mbp","Plp1","Cxcrl1","Inpp5d"))
dev.off()
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
