#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/MeA_Inhibexcitheat2_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)



mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")

pdf(file=paste0(output,"labeledUMAPFigS4.pdf"))
DimPlot(mySeurat, label=TRUE, reduction="umap",label.size=6)
dev.off()
pdf(file=paste0(output,"Gad1.pdf"))
VlnPlot(mySeurat, features=c("Gad1"), idents=c("5","22","31"), pt.size=0,  cols=c("gray","gray","gray"))
dev.off()
pdf(file=paste0(output,"vglut2.pdf"))
VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("5","22","31"), pt.size=0,  cols=c("gray","gray","gray"))
dev.off()
gad1f <- VlnPlot(mySeurat, features=c("Gad1"), idents=c("5","22","31"), pt.size=0,  cols=c("gray","gray","gray"))
gad1fdat <- gad1f$data
write.csv(gad1fdat, file=paste0(output, "gad1mixedvlndata.csv"))
pdf(file=paste0(output,"mixedgad1cust.pdf"))
ggplot(gad1fdat, aes(x=Gad1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()

ebf3 <- VlnPlot(mySeurat, features=c("Ebf3"), idents=c("22"), pt.size=0,  cols=c("gray","gray","gray"))
ebf3 <- ebf3$data
write.csv(ebf3, file=paste0(output, "ebf3mixedvlndata.csv"))
pdf(file=paste0(output,"mixedebf3cust.pdf"))
ggplot(ebf3, aes(x=Ebf3, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=FALSE, fill="black", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()


glutf <- VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("5","22","31"), pt.size=0,  cols=c("gray","gray","gray"))
glutfdat <- glutf$data
write.csv(glutfdat, file=paste0(output, "glutmixedvlndata.csv"))
pdf(file=paste0(output,"mixedglutcust.pdf"))
ggplot(glutfdat, aes(x=Slc17a6, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=FALSE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()


gad1f <- VlnPlot(mySeurat, features=c("Gad1"), idents=c("21"), pt.size=0,  cols=c("gray","gray","gray"))
gad1fdat <- gad1f$data
write.csv(gad1fdat, file=paste0(output, "gad1femalevlndata.csv"))
pdf(file=paste0(output,"femalegad1cust.pdf"))
ggplot(gad1fdat, aes(x=Gad1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()

glutf <- VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("21"), pt.size=0,  cols=c("gray","gray","gray"))
glutfdat <- glutf$data
write.csv(glutfdat, file=paste0(output, "glutfemalevlndata.csv"))
pdf(file=paste0(output,"femaleglutcust.pdf"))
ggplot(glutfdat, aes(x=Slc17a6, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()


pdf(file=paste0(output,"mixedGad1.pdf"))
VlnPlot(mySeurat, features=c("Gad1"), idents=c("5","23","32"), pt.size=0, cols=c("gray","gray","gray"))
dev.off()
pdf(file=paste0(output,"femaleupGad1.pdf"))
VlnPlot(mySeurat, features=c("Gad1"), idents=c("21"), pt.size=0, cols=c("gray"))
dev.off()
pdf(file=paste0(output,"femaleupvmat.pdf"))
VlnPlot(mySeurat, features=c("Slc18a2"), idents=c("21"), pt.size=0, cols=c("gray"))
dev.off()
pdf(file=paste0(output,"mixedUpSlc17a6.pdf"))
VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("5","23","32"), pt.size=0, cols=c("gray","gray","gray"))
dev.off()
pdf(file=paste0(output,"femaleupSlc17a6.pdf"))
VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("21"), pt.size=0, cols=c("gray"))
dev.off()
pdf(file=paste0(output,"inhibexcit.pdf"))
DimPlot(mySeurat, cols=c("0"="blue", "1"="red", "2"="red", "3"="blue", "4"="green", "5"="red", "6"="red", 
	"7"="red", "8"="blue", "9"="green", "10"="red", "11"="red", "12"="blue", "13"="blue", 
	"14"="blue", "15"="red", "16"="blue", "17"="blue", "18"="red", "19"="red", "20"="blue", 
	"21"="blue", "22"="blue", "23"="red", "24"="red", "25"="green", "26"="red", "27"="red", "28"="blue", "29"="red", "30"="blue", "31"="blue", "32"="blue", "33"="blue", "34"="red"), label=FALSE, reduction="umap")
dev.off() 
genelist <- read.table("Genelists/MeA_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
MvP <- read.table("Genelists/MeA_MvP_1.5cutoff.txt", header=FALSE)
MvU <- read.table("Genelists/MeA_MvU_1.5cutoff.txt", header=FALSE)
sDEGs <- rbind(MvP, MvU)
sDEGs <- unlist(sDEGs)
sDEGs <- unique(sDEGs)
filtered.sDEGs <- intersect(sDEGs, genes.10x)
all <- rbind(genelist, MvP, MvU)
all <- unlist(all)
all <- unique(all)
filtered.all <- intersect(all, genes.10x)
data <- mySeurat@assays[["RNA"]]@data

AvgAllDEGs <- AverageExpression(mySeurat, features=filtered.all, assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/AllDEGsavg.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/AllDEGsavg.txt", header=TRUE, row.names=1)
data <- as.matrix(data)
dim(data)
head(data)
data <- scale(t(data))
idents <- mySeurat@active.ident
idents <- as.matrix(idents)
head(idents)
pdf(file=paste0(output,"VMH_sDEGheatmaptrowclustonly.pdf"), width=30, height=30)
pheatmap(data, cluster_cols=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"VMH_sDEGheatmapt.pdf"), width=30, height=30)
pheatmap(data, border_color=NA, color=viridis(500))
dev.off()
idents <- mySeurat@active.ident
df <- cbind(data, idents)
dim(df)
head(df)
DEGs <- df[,filtered.genelist]
dim(DEGs)
head(DEGs)
mySeurat <- ScaleData(mySeurat, features=rownames(mySeurat))
pdf(file=paste0(output,"PCAClust_sDEGheatmap.pdf"))
DoHeatmap(mySeurat, features=filtered.genelist)
dev.off()
