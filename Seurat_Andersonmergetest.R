#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/Anderson/AndersonFullMerge_removenoisy"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Rep1.data <- Read10X(data.dir = "2018_0410")
GroupRep1 <- CreateSeuratObject(counts = Rep1.data, project = "GroupRep1", min.cells=3, min.features=200)
Rep2.data <- Read10X(data.dir= "2018_0420")
GroupRep2 <- CreateSeuratObject(counts = Rep2.data, project = "GroupRep2", min.cells=3, min.features=200)
Rep3.data <- Read10X(data.dir= "2018_0429")
GroupRep3 <- CreateSeuratObject(counts = Rep3.data, project = "GroupRep3", min.cells=3, min.features=200)
Rep4.data <- Read10X(data.dir= "2018_0520")
GroupRep4 <- CreateSeuratObject(counts = Rep4.data, project = "GroupRep4", min.cells=3, min.features=200)
Rep5.data <- Read10X(data.dir= "2018_0527")
GroupRep5 <- CreateSeuratObject(counts = Rep5.data, project = "GroupRep5", min.cells=3, min.features=200)
RepS1.data <- Read10X(data.dir = "2019_0330")
SingleRep1 <- CreateSeuratObject(counts = RepS1.data, project = "SingleRep1", min.cells=3, min.features=200)
RepS2.data <- Read10X(data.dir= "2018_0812_1")
SingleRep2 <- CreateSeuratObject(counts = RepS2.data, project = "SingleRep2", min.cells=3, min.features=200)
RepG1.data <- Read10X(data.dir= "2019_0325")
RepG1 <- CreateSeuratObject(counts = RepG1.data, project = "GroupFear", min.cells=3, min.features=200)
GroupRep1$Condition <- "Control"
GroupRep2$Condition <- "Control"
GroupRep3$Condition <- "Control"
GroupRep4$Condition <- "Control"
GroupRep5$Condition <- "Control"
SingleRep1$Condition <- "Single Housed Fear"
SingleRep2$Condition <- "Single Housed Fear"
RepG1$Condition <- "Group Housed Fear"
mySeurat <- merge(GroupRep1, y=c(GroupRep2, GroupRep3, GroupRep4, GroupRep5, SingleRep1, SingleRep2, RepG1), add.cell.ids=c("GroupRep1","GroupRep2", "GroupRep3","GroupRep4","GroupRep5","SingleRep1","SingleRep2","RepG1"), project="FullHouseMerge")
mySeurat
mySeurat[["percent.mt"]] <- PercentageFeatureSet(mySeurat, pattern = "^mt-")
mySeurat <- subset(mySeurat, subset = nFeature_RNA >600 & nCount_RNA < 30000 & percent.mt < 15)
head(mySeurat[[]])
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/AndersonExclude.txt", header=FALSE)
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
mySeurat <- ScaleData(mySeurat, vars.to.regress="percent.mt")
mySeurat <- RunPCA(mySeurat, features=hvg.final)
pdf(paste0(output,"_PCAplot.pdf"))
DimPlot(mySeurat, reduction="pca", group.by="orig.ident")
dev.off()
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
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
pdf(paste0(output,"_UMAP_condsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Condition")
dev.off()
pdf(paste0(output,"_UMAP_condlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Condition")
dev.off()
mySeurat <- RunTSNE(mySeurat, features=VariableFeatures(object=mySeurat))
pdf(paste0(output,"_tsne.pdf"))
DimPlot(mySeurat, reduction="tsne", label=TRUE)
dev.off()
pdf(paste0(output,"_tsne_origsplit.pdf"))
DimPlot(mySeurat, reduction="tsne", label=TRUE, split.by="orig.ident")
dev.off()
pdf(paste0(output,"_tsne_origlabel.pdf"))
DimPlot(mySeurat, reduction="tsne", label=TRUE, group.by="orig.ident")
dev.off()
pdf(paste0(output,"_tsne_condsplit.pdf"))
DimPlot(mySeurat, reduction="tsne", label=TRUE, split.by="Condition")
dev.off()
pdf(paste0(output,"_tsne_condlabel.pdf"))
DimPlot(mySeurat, reduction="tsne", label=TRUE, group.by="Condition")
dev.off()

