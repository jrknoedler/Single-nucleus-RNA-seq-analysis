#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_Anderson/AndersonFemale_GliaExcitPgrFilt"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)


mySeurat <- readRDS("Seurat/VMH_Anderson/AndersonFemale_GliaFilt7.rds")
mySeurat <- subset(mySeurat, idents=c("0","3","11","17","18","21","33","39"))
#mySeurat[["percent.mt"]] <- PercentageFeatureSet(mySeurat, pattern = "mt-^")

mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.mt"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
IEGs <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/AndersonExclude.txt", header=FALSE)
IEGs <- unlist(IEGs)
genelist
genelist <- unlist(genelist)
genelist
exclude <- rbind(genelist,IEGs)
exclude <- as.matrix(exclude)

mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, exclude)
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
VlnPlot(mySeurat, features=c("Slc32a1","Slc17a6","Esr1","Cckar"), ncol=1, pt.size=0)
dev.off()
SeuratScaled <- ScaleData(mySeurat)
pdf(paste0(output,"_Markerdot.pdf"), width=40, height=40)
DotPlot(SeuratScaled, features=c("Gad1","Slc17a6","Esr1","Cckar","Mbp","Plp1","Cxcr1l","Inpp5d","Iqschfp"))
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
