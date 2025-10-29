#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"


library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)
library(patchwork)
library(Matrix)

BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
pdf(file="Seurat/ReclusteredBNST_Exciteinhib.pdf", width=20)
VlnPlot(BNST, features=c("Gad1","Slc17a6"), pt.size=0, ncol=1)
dev.off()

Matrix::colMedians(mySeurat, slot="counts")

DefaultAssay(BNST) <- "RNA"
BNST <- NormalizeData(BNST)
BNST

BNSTgenes <- read.table("topGO/Total_ByRegion/BNST_Genesonly.txt", header=FALSE)

BNSTgenes <- unlist(BNSTgenes)

genes.BNST <- (x=rownames(x=BNST))
unlist(genes.BNST)
filtered.BNST <- intersect(BNSTgenes, genes.BNST)
filtered.BNST


AvgAllDEGs <- AverageExpression(BNST, features=filtered.BNST, assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/DEGheatmapfinal.csv")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/DEGheatmapfinal.csv", header=TRUE, row.names=1)
data <- as.matrix(data)

data <- scale(t(data))

BNSTmap <- pheatmap(data, cluster_rows=FALSE)

data <- data[,BNSTmap$tree_col[["order"]]]
write.csv(data, file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/ReorderedDEGheatmapfinal.csv")

VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
pdf(file="Seurat/ReclusteredVMH_Exciteinhib.pdf", width=20)
VlnPlot(VMH, features=c("Gad1","Slc17a6"), pt.size=0, ncol=1)
dev.off()
DefaultAssay(VMH) <- "RNA"
VMH <- NormalizeData(VMH)
VMH
VMHgenes <- read.table("topGO/Total_ByRegion/VMH_Genesonly.txt", header=FALSE)


VMHgenes <- unlist(VMHgenes)

genes.VMH <- (x=rownames(x=VMH))
unlist(genes.VMH)
filtered.VMH <- intersect(VMHgenes, genes.VMH)
filtered.VMH

AvgAllDEGs <- AverageExpression(VMH, features=filtered.VMH, assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/DEGheatmapfinal.csv")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/DEGheatmapfinal.csv", header=TRUE, row.names=1)
data <- as.matrix(data)

data <- scale(t(data))

VMHmap <- pheatmap(data, cluster_rows=FALSE)

data <- data[,VMHmap$tree_col[["order"]]]
write.csv(data, file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/ReorderedDEGheatmapfinal.csv")


POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
pdf(file="Seurat/ReclusteredPOA_Exciteinhib.pdf", width=20)
VlnPlot(POA, features=c("Gad1","Slc17a6"), pt.size=0, ncol=1)
dev.off()
POA
POAgenes <- read.table("topGO/Total_ByRegion/POA_Genesonly.txt", header=FALSE)
DefaultAssay(POA) <- "RNA"
POA <- NormalizeData(POA)

POAgenes <- unlist(POAgenes)

genes.POA <- (x=rownames(x=POA))
unlist(genes.POA)
filtered.POA <- intersect(POAgenes, genes.POA)
filtered.POA

AvgAllDEGs <- AverageExpression(POA, features=filtered.POA, assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/DEGheatmapfinal.csv")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/POA/DEGheatmapfinal.csv", header=TRUE, row.names=1)
data <- as.matrix(data)

data <- scale(t(data))
POAmap <- pheatmap(data, cluster_rows=FALSE)

data <- data[,POAmap$tree_col[["order"]]]
write.csv(data, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/ReorderedDEGheatmapfinal.csv")




MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")

MeA

MeAgenes <- read.table("topGO/Total_ByRegion/MeA_Genesonly.txt", header=FALSE)

MeAgenes <- unlist(MeAgenes)

genes.MeA <- (x=rownames(x=MeA))
unlist(genes.MeA)
filtered.MeA <- intersect(MeAgenes, genes.MeA)
filtered.MeA
DefaultAssay(MeA) <- "RNA"
MeA <- NormalizeData(MeA)
pdf(file="Seurat/ReclusteredMeA_Exciteinhib.pdf", width=20)
VlnPlot(MeA, features=c("Gad1","Slc17a6"), pt.size=0, ncol=1)
dev.off()

AvgAllDEGs <- AverageExpression(MeA, features=filtered.MeA, assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/DEGheatmapfinal.csv")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/DEGheatmapfinal.csv", header=TRUE, row.names=1)
data <- as.matrix(data)

data <- scale(t(data))
MeAmap <- pheatmap(data, cluster_rows=FALSE)

data <- data[,MeAmap$tree_col[["order"]]]
write.csv(data, file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/ReorderedDEGheatmapfinal.csv")
