#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_Merge_Inhibexcitheat_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)



mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")

pdf(file=paste0(output,"maleupGad1.pdf"))
VlnPlot(mySeurat, features=c("Gad1"), idents=c("19"), pt.size=0, cols=c("gray"))
dev.off()
gad <- VlnPlot(mySeurat, features=c("Gad1"), idents=c("19"), pt.size=0)
gaddat <- gad$data
write.csv(gaddat, file=paste0(output, "gadvlndata.csv"))
pdf(file=paste0(output,"malegad1cust.pdf"))
ggplot(gaddat, aes(x=Gad1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,3)
dev.off()
pdf(file=paste0(output,"maleupSlc17a6.pdf"))
VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("19"), pt.size=0, cols=c("gray"))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("19"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "glustvlndata.csv"))
pdf(file=paste0(output,"maleglutcust.pdf"))
ggplot(glutdat, aes(x=Slc17a6, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()
pdf(file=paste0(output,"maleupCyp19a1.pdf"))
VlnPlot(mySeurat, features=c("Cyp19a1"), idents=c("19"), pt.size=0, cols=c("gray"))
dev.off()
cyp <- VlnPlot(mySeurat, features=c("Cyp19a1"), idents=c("19"), pt.size=0)
cypdat <- cyp$data
write.csv(cypdat, file=paste0(output, "glustvlndata.csv"))
pdf(file=paste0(output,"malecypcust.pdf"))
ggplot(cypdat, aes(x=Cyp19a1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()


gad1f <- VlnPlot(mySeurat, features=c("Gad1"), idents=c("8","21","26"), pt.size=0,  cols=c("gray","gray","gray"))
gad1fdat <- gad1f$data
write.csv(gad1fdat, file=paste0(output, "gad1fvlndata.csv"))
pdf(file=paste0(output,"femalegad1cust.pdf"))
ggplot(gad1fdat, aes(x=Gad1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()

glutf <- VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("8","21","26"), pt.size=0,  cols=c("gray","gray","gray"))
glutfdat <- glutf$data
write.csv(glutfdat, file=paste0(output, "glutfvlndata.csv"))
pdf(file=paste0(output,"femaleglutcust.pdf"))
ggplot(glutfdat, aes(x=Slc17a6, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()


vmatf <- VlnPlot(mySeurat, features=c("Th"), idents=c("8","21","26"), pt.size=0,  cols=c("gray","gray","gray"))
vmatfdat <- vmatf$data
write.csv(vmatfdat, file=paste0(output, "thfvlndata.csv"))
pdf(file=paste0(output,"femaleThcust.pdf"))
ggplot(vmatfdat, aes(x=Th, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,1)
dev.off()


pdf(file=paste0(output,"labeledUMAPFigS4.pdf"))
DimPlot(mySeurat, label=TRUE, reduction="umap",label.size=6)
dev.off()
pdf(file=paste0(output,"sexbias.pdf"))
DimPlot(mySeurat, cols=c("0"="red", "1"="blue", "2"="blue", "3"="blue", "4"="blue", "5"="red", "6"="blue", 
	"7"="red", "8"="blue", "9"="blue", "10"="blue", "11"="red", "12"="red", "13"="red", 
	"14"="blue", "15"="red", "16"="blue", "17"="blue", "18"="red", "19"="blue", "20"="blue", 
	"21"="blue", "22"="red", "23"="blue", "24"="blue", "25"="blue", "26"="blue", "27"="red"), label=FALSE, reduction="umap")
dev.off() 
genelist <- read.table("Genelists/VMH_PvU_1.5.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
MvP <- read.table("Genelists/VMH_MvP_1.5.txt", header=FALSE)
MvU <- read.table("Genelists/VMH_MvU_1.5.txt", header=FALSE)
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
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/AllDEGsavg.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/AllDEGsavg.txt", header=TRUE, row.names=1)
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

