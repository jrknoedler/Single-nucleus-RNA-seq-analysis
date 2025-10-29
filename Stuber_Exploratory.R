#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/Stuber_Habenula/Stuber_Round1"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)

Control.data <- Read10X(data.dir="Stuber_Habenula/control")
Control <- CreateSeuratObject(counts = Control.data, project = "Control", min.cells=3, min.features=200)
Control$condition <- "Control"

Footshock.data <- Read10X(data.dir="Stuber_Habenula/footshock")
Footshock <- CreateSeuratObject(counts = Footshock.data, project = "Footshock", min.cells=3, min.features=200)
Footshock$condition <- "Footshock"

mySeurat <- merge(Control, y=c(Footshock), add.cell.ids=c("Control","Footshock"), project="StuberMerged")
mySeurat
mySeurat[["percent.mito"]] <- PercentageFeatureSet(mySeurat, pattern = "mt-^")
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident", "percent.mito"))
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
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_condlabel.pdf"))
DimPlot(mySeurat, reduction="umap", group.by="condition", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_condsplit.pdf"))
DimPlot(mySeurat, reduction="umap", split.by="condition", label=TRUE)
dev.off()
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=60)
VlnPlot(mySeurat, features=c("Slc32a1","Slc17a6","Tacr1","Mbp","Plp1","Cxcr1","Inppd5"), ncol=1, pt.size=0)
dev.off()
SeuratScaled <- ScaleData(mySeurat)
pdf(paste0(output,"_Markerdot.pdf"), width=40, height=40)
DotPlot(SeuratScaled, features=c("Slc32a1","Slc17a6","Tacr1","Mbp","Plp1","Cxcr1","Inppd5"))
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=paste0(output, ".rds"))