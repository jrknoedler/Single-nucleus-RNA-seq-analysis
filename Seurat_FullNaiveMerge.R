#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/WinterBreakMergetest/FullNaiveMerge"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
MaleBNST.data <- Read10X(data.dir = "MaleBNST_Runscombined/outs/filtered_feature_bc_matrix")
MaleBNST <- CreateSeuratObject(counts = MaleBNST.data, project = "MaleBNST", min.cells=3, min.features=200)
MalePOA.data <- Read10X(data.dir= "MalePOAIntact/outs/filtered_feature_bc_matrix")
MalePOA <- CreateSeuratObject(counts = MalePOA.data, project = "MalePOA", min.cells=3, min.features=200)
PrimedBNST.data <- Read10X(data.dir= "MaleIntactVMH_retry/outs/filtered_feature_bc_matrix")
PrimedBNST <- CreateSeuratObject(counts = PrimedBNST.data, project = "PrimedBNST", min.cells=3, min.features=200)
PrimedPOA.data <- Read10X(data.dir= "MaleMeAIntact/outs/filtered_feature_bc_matrix")
PrimedPOA <- CreateSeuratObject(counts = PrimedPOA.data, project = "PrimedPOA", min.cells=3, min.features=200)
MaleMeA.data <- Read10X(data.dir= "FemaleBNSTPrimed_Putative_retry/outs/filtered_feature_bc_matrix")
MaleMeA <- CreateSeuratObject(counts = MaleMeA.data, project = "MaleMeA", min.cells=3, min.features=200)
MaleVMH.data <- Read10X(data.dir = "FemalePOAPrimed_Putative2/outs/filtered_feature_bc_matrix")
MaleVMH <- CreateSeuratObject(counts = MaleVMH.data, project = "MaleVMH", min.cells=3, min.features=200)
PrimedMeA.data <- Read10X(data.dir= "FemaleMeAPrimed/outs/filtered_feature_bc_matrix")
PrimedMeA <- CreateSeuratObject(counts = PrimedMeA.data, project = "PrimedMeA", min.cells=3, min.features=200)
PrimedVMH.data <- Read10X(data.dir= "FemaleVMHPrimed/outs/filtered_feature_bc_matrix")
PrimedVMH <- CreateSeuratObject(counts = PrimedVMH.data, project = "PrimedVMH", min.cells=3, min.features=200)
UnprimedBNST.data <- Read10X(data.dir= "FemaleBNSTUnprimed/outs/filtered_feature_bc_matrix")
UnprimedBNST <- CreateSeuratObject(counts = UnprimedBNST.data, project = "UnprimedBNST", min.cells=3, min.features=200)
UnprimedPOA.data <- Read10X(data.dir= "FemalePOAUnprimed/outs/filtered_feature_bc_matrix")
UnprimedPOA <- CreateSeuratObject(counts = UnprimedPOA.data, project = "UnprimedPOA", min.cells=3, min.features=200)
MaleVMH$Sex <- "Male"
MaleVMH$Region <- "VMH"
MaleVMH$Batch <- "3"
MaleVMH$Hormone <- "Intact"
MaleBNST$Sex <- "Male"
MaleBNST$Region <- "BNST"
MaleBNST$Batch <- "1"
MaleBNST$Hormone <- "Intact"
MalePOA$Sex <- "Male"
MalePOA$Region <- "POA"
MalePOA$Batch <- "2"
MalePOA$Hormone <- "Intact"
MaleMeA$Sex <- "Male"
MaleMeA$Region <- "MeA"
MaleMeA$Batch <- "3"
MaleMeA$Hormone <- "Intact"
PrimedBNST$Sex <- "Female"
PrimedBNST$Region <- "BNST"
PrimedBNST$Batch <- "3"
PrimedBNST$Hormone <- "Primed"
PrimedPOA$Sex <- "Female"
PrimedPOA$Region <- "POA"
PrimedPOA$Batch <- "3"
PrimedPOA$Hormone <- "Primed"
PrimedMeA$Sex <- "Female"
PrimedMeA$Region <- "MeA"
PrimedMeA$Batch <- "3"
PrimedMeA$Hormone <- "Primed"
PrimedVMH$Sex <- "Female"
PrimedVMH$Region <- "VMH"
PrimedVMH$Batch <- "3"
PrimedVMH$Hormone <- "Primed"
UnprimedPOA$Sex <- "Female"
UnprimedPOA$Region <- "POA"
UnprimedPOA$Batch <- "2"
UnprimedPOA$Hormone <- "Unprimed"
UnprimedBNST$Sex <- "Female"
UnprimedBNST$Region <- "BNST"
UnprimedBNST$Batch <- "2"
UnprimedBNST$Hormone <- "Unprimed"
mySeurat <- merge(MaleBNST, y=c(MalePOA, MaleMeA, MaleVMH, UnprimedPOA, UnprimedBNST, PrimedPOA, PrimedBNST, PrimedMeA, PrimedVMH), add.cell.ids=c("MaleBNST","MalePOA", "MaleMeA","MaleVMH","PrimedBNST","PrimedPOA","PrimedMeA","PrimedVMH","UnprimedBNST","UnprimedPOA"), project="FullHouseMerge")
mySeurat
mySeurat[["percent.mt"]] <- PercentageFeatureSet(mySeurat, pattern = "^mt-")
mySeurat <- subset(mySeurat, subset = nFeature_RNA >600 & nCount_RNA < 60000 & percent.mt < 15)
head(mySeurat[[]])
mySeurat <- NormalizeData(mySeurat)
mySeurat <- FindVariableFeatures(mySeurat, selection.method="vst", nfeatures=2000)
mySeurat <- ScaleData(mySeurat)
mySeurat <- RunPCA(mySeurat, npcs=50)
pdf(paste0(output,"_PCAplot.pdf"))
DimPlot(mySeurat, reduction="pca", group.by="orig.ident")
dev.off()
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:50)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:50)
mySeurat <- FindClusters(mySeurat, resolution=1)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_origsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="orig.ident")
dev.off()
pdf(paste0(output,"_UMAP_origlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="orig.ident")
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_regionsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Region")
dev.off()
pdf(paste0(output,"_UMAP_regionlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Region")
dev.off()
pdf(paste0(output,"_UMAP_batchsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Batch")
dev.off()
pdf(paste0(output,"_UMAP_batchlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Batch")
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
FeaturePlot(mySeurat, features=c("Xist"), split.by="orig.ident")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
