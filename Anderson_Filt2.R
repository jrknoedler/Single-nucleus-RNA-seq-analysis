#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_Anderson/Anderson_NonneuronalFilt2"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)


mySeurat <- readRDS("Seurat/VMH_Anderson/Anderson_GliaFilt.rds")
mySeurat <- subset(mySeurat, idents=c("17","18","22","29","34","35","36","44","55","56","57","59","60"), invert=TRUE)
mySeurat[["percent.mt"]] <- PercentageFeatureSet(mySeurat, pattern = "mt-^")

mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.mt"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
#genelist <- as.matrix(genelist)
genelist
genelist2 <- read.table("Andersonexclude.txt", header=FALSE)
genelist2 <- unlist(genelist2)
genelistcomb <- rbind(genelist, genelist2)
genelistcomb <- unique(genelistcomb)
genelistcomb <- as.matrix(genelistcomb)
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelistcomb)
hvg.final
mySeurat <- RunPCA(mySeurat, features=hvg.final)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
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
DotPlot(SeuratScaled, features=c("Gad1","Slc17a6","Esr1","Cckar","Mbp","Plp1","Cxcr1","Inppd5"))
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
