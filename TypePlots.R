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

BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
pdf(file="Seurat/ReclusteredBNST_Exciteinhib.pdf", width=20)
VlnPlot(BNST, features=c("Gad1","Slc17a6"), pt.size=0, ncol=1)
dev.off()
DefaultAssay(BNST) <- "RNA"
BNST <- NormalizeData(BNST)
BNST
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(BNST)
BNST <- RenameIdents(BNST, new.cluster.ids)

BNSTgenes <- read.table("topGO/Total_ByRegion/BNST_Genesonly.txt", header=FALSE)
p <- DimPlot(BNST, reduction="umap", label=TRUE, label.size=8)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
p2 <- DimPlot(BNST, reduction="umap", label=FALSE, cols=c("1"="red","2"="red","3"="red","4"="red","5"="red","6"="light blue","7"="red","8"="red","9"="red","10"="red","11"="red","12"="red","13"="red","14"="red","15"="light blue","16"="red","17"="red","18"="red","19"="red","20"="light blue","21"="red","22"="red","23"="red","24"="red","25"="red","26"="red","27"="red","28"="red","29"="red","30"="red","31"="red","32"="red","33"="red","34"="red","35"="light blue","36"="red"))
p2f <- p2 +theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
pdf(file="Seurat/BNSTtypeplot.pdf", width=12)
pd+p2f
dev.off()
counts <- table(Idents(BNST), BNST$Hormone)
write.csv(counts, file=paste0("Seurat/BNST_finalcountsperclust.csv"))


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
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(VMH)
VMH <- RenameIdents(VMH, new.cluster.ids)

p <- DimPlot(VMH, reduction="umap", label=TRUE, label.size=8)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
p2 <- DimPlot(VMH, reduction="umap", label=FALSE, cols=c("1"="red","2"="light blue","3"="light blue","4"="light blue","5"="red","6"="light blue","7"="light blue","8"="red","9"="red","10"="red","11"="light blue","12"="light blue","13"="red","14"="red","15"="light blue","16"="light blue","17"="light blue","18"="red","19"="light blue","20"="light blue","21"="light blue","22"="red","23"="red","24"="light blue","25"="light blue","26"="light blue","27"="light blue"))
p2f <- p2 +theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
pdf(file="Seurat/VMHtypeplot.pdf", width=12)
pd+p2f
dev.off()
counts <- table(Idents(VMH), VMH$Hormone)
write.csv(counts, file=paste0("Seurat/VMH_finalcountsperclust.csv"))

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
pdf(file="Seurat/ReclusteredPOA_Exciteinhib.pdf", width=12)
VlnPlot(POA, features=c("Gad1","Slc17a6"), pt.size=0, ncol=1)
dev.off()
POA
POAgenes <- read.table("topGO/Total_ByRegion/POA_Genesonly.txt", header=FALSE)
DefaultAssay(POA) <- "RNA"
POA <- NormalizeData(POA)

POAgenes <- unlist(POAgenes)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(POA)
POA <- RenameIdents(POA, new.cluster.ids)

p <- DimPlot(POA, reduction="umap", label=TRUE, label.size=8)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
p2 <- DimPlot(POA, reduction="umap", label=FALSE, cols=c("1"="light blue","2"="red","3"="red","4"="light blue","5"="red","6"="light blue","7"="red","8"="red","9"="red","10"="light blue","11"="red","12"="light blue","13"="red","14"="red","15"="red","16"="light blue","17"="red","18"="light blue","19"="red","20"="light blue","21"="red","22"="red","23"="light blue","24"="red","25"="light blue","26"="red","27"="light blue","28"="red","29"="red","30"="red","31"="red","32"="red","33"="light blue","34"="red","35"="red","36"="red", "37"="red","38"="red","39"="light blue","40"="light blue"))
p2f <- p2 +theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
pdf(file="Seurat/POAtypeplot.pdf", width=12)
pd+p2f
dev.off()
counts <- table(Idents(POA), POA$Hormone)
write.csv(counts, file=paste0("Seurat/POA_finalcountsperclust.csv"))


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

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
names(new.cluster.ids) <- levels(MeA)
MeA <- RenameIdents(MeA, new.cluster.ids)

p <- DimPlot(MeA, reduction="umap", label=TRUE, label.size=8)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
p2 <- DimPlot(MeA, reduction="umap", label=FALSE, cols=c("1"="light blue","2"="red","3"="red","4"="light blue","5"="red","6"="red","7"="light green","8"="red","9"="light blue","10"="light green","11"="red","12"="red","13"="light blue","14"="light blue","15"="light green","16"="light blue","17"="red","18"="light blue","19"="red","20"="red","21"="light blue","22"="light blue","23"="light blue","24"="light blue","25"="red","26"="red","27"="red","28"="red","29"="light blue","30"="light blue","31"="light blue","32"="light blue","33"="light blue","34"="red"))
p2f <- p2 +theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
pdf(file="Seurat/MeAtypeplot.pdf", width=12)
pd+p2f
dev.off()
counts <- table(Idents(MeA), MeA$Hormone)
write.csv(counts, file=paste0("Seurat/MeA_finalcountsperclust.csv"))