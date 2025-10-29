#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_Main6_altSCTreclust_linevln"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_SCTFiltered1_20pcs_0.8.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mylevels <- c(6,10,20,33,0,1,2,3,4,5,7,8,9,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32)

pdf(file=paste0(output,"Tacr1_peptides.pdf"), height=20, width=25)
VlnPlot(mySeurat, features=c("Tacr1","Cartpt","Cck","Sst","Vip","Tac2"), ncol=1, pt.size=0)
dev.off()

glut <- VlnPlot(mySeurat, features=c("Esr1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Esr1vlndata.csv"))
pdf(file=paste0(output,"Esr1vln.pdf"), width=40)
ggplot(glutdat, aes(y=Esr1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", scale="width") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank())
dev.off()

glut <- VlnPlot(mySeurat, features=c("Gad1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Gad1vlndata.csv"))
pdf(file=paste0(output,"Gad1vln.pdf"), width=40)
ggplot(glutdat, aes(y=Gad1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", scale="width") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()

glut <- VlnPlot(mySeurat, features=c("Cyp19a1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Cyp19a1vlndata.csv"))
pdf(file=paste0(output,"Cyp19a1vln.pdf"), width=40)
ggplot(glutdat, aes(y=Cyp19a1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", scale="width") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()

glut <- VlnPlot(mySeurat, features=c("Slc17a6"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Slc17a6vlndata.csv"))
pdf(file=paste0(output,"Slc17a6vln.pdf"), width=30)
ggplot(glutdat, aes(y=Slc17a6, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", scale="width") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()

glut <- VlnPlot(mySeurat, features=c("Tac1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Tac1vlndata.csv"))
pdf(file=paste0(output,"Tac1vln.pdf"), width=30)
ggplot(glutdat, aes(y=Tac1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="gray", scale="width") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) 
dev.off()
SFARI <- read.table("Heatmaps/NewSFARI.txt")
SFARI <- unlist(SFARI)
SFARI <- as.matrix(SFARI)
AvgSFARI <- AverageExpression(mySeurat, features=SFARI, assays=c("RNA"))
AvgAllDEGs <- AverageExpression(mySeurat, features=c("Tac1","Cyp19a1","Esr2","Sox5","Moxd1","Pappa","Greb1","Cck","Cartpt"), assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/Tac1markerssavg2.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/Tac1markerssavg2.txt", header=TRUE, row.names=1)
data <- as.matrix(data)
dim(data)
head(data)
data <- scale(t(data))
write.csv(data, file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/Tac1_scaledredo.csv")
idents <- mySeurat@active.ident
idents <- as.matrix(idents)
head(idents)
pdf(file=paste0(output,"Tac1heatmaptrowclustonly.pdf"))
pheatmap(data, cluster_cols=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"Tac1heatmapnoclust.pdf"))
pheatmap(data, cluster_cols=FALSE, cluster_rows=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"Tac1heatmapt.pdf"), width=30, height=30)
pheatmap(data, border_color=NA, color=viridis(500))
dev.off()
write.csv(AvgSFARI, file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/AllSFARI.csv")
data2 <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/AllSFARI.csv", header=TRUE, row.names=1)
data2 <- as.matrix(data2)
data2 <- scale(t(data2))
data2
pdf(file=paste0(output,"SFARIsDEGs.pdf"))
pheatmap(data2, color=viridis(500))
dev.off()
pdf(file=paste0(output,"SFARIsDEGsDEGclust.pdf"))
pheatmap(data2, cluster_rows=FALSE, color=viridis(500))
dev.off()