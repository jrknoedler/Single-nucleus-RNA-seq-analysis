#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_Inhibexcitheat_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)



mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)

cyp <- VlnPlot(mySeurat, features=c("Cyp19a1"), idents=c("10","13","19","24","33"), pt.size=0)
cypdat <- cyp$data
write.csv(cypdat, file=paste0(output, "glustvlndata.csv"))
pdf(file=paste0(output,"malecypcust.pdf"))
ggplot(cypdat, aes(x=Cyp19a1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,3)
dev.off()
pdf(file=paste0(output,"Th.pdf"))
VlnPlot(mySeurat, features=c("Th"), idents=c("21","26"), pt.size=0,  cols=c("gray","gray","gray"))
dev.off()
vmatf <- VlnPlot(mySeurat, features=c("Th"), idents=c("21","26"), pt.size=0,  cols=c("gray","gray","gray"))
vmatfdat <- vmatf$data
write.csv(vmatfdat, file=paste0(output, "thfvlndata.csv"))
pdf(file=paste0(output,"femalethcust.pdf"))
ggplot(vmatfdat, aes(x=Th, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,1)
dev.off()


pdf(file=paste0(output,"femaleupGad1.pdf"))
VlnPlot(mySeurat, features=c("Gad1"), idents=c("21","26"), pt.size=0, cols=c("gray","gray"))
dev.off()

gadf <- VlnPlot(mySeurat, features=c("Gad1"), idents=c("21","26"), pt.size=0, cols=c("gray","gray"))
gadfdat <- gadf$data
pdf(file=paste0(output,"femalegad1cust.pdf"))
ggplot(gadfdat, aes(x=Gad1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()
write.csv(gadfdat, file=paste0(output, "gadfvlndata.csv"))
pdf(file=paste0(output,"maleupGad1.pdf"))
VlnPlot(mySeurat, features=c("Gad1"), idents=c("10","13","19","24","33"), pt.size=0, cols=c("gray","gray","gray","gray","gray"))
dev.off()
gad <- VlnPlot(mySeurat, features=c("Gad1"), idents=c("10","13","19","24","33"), pt.size=0, cols=c("gray","gray","gray","gray","gray"))
gaddat <- gad$data
write.csv(gaddat, file=paste0(output, "gadvlndata.csv"))
pdf(file=paste0(output,"malecat1cust.pdf"))
ggplot(gaddat, aes(x=Gad1, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,3)
dev.off()
pdf(file=paste0(output,"maleupCyp19a1.pdf"))
VlnPlot(mySeurat, features=c("Cyp19a1"), idents=c("10","13","19","24","33"), pt.size=0, cols=c("gray","gray","gray","gray","gray"))
dev.off()
pdf(file=paste0(output,"maleupSlc17a6.pdf"))
VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("10","13","19","24","33"), pt.size=0, cols=c("gray","gray","gray","gray","gray"))
dev.off()
glut <- VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("10","13","19","24","33"), pt.size=0, cols=c("gray","gray","gray","gray","gray"))
glutdat <- glut$data
pdf(file=paste0(output,"maleglutcust.pdf"))
ggplot(glutdat, aes(x=Slc17a6, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2)
dev.off()
write.csv(glutdat, file=paste0(output, "glutvlndata.csv"))
pdf(file=paste0(output,"femaleupSlc17a6.pdf"))
VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("21","26"), pt.size=0, cols=c("gray","gray"))
dev.off()
glutf <- VlnPlot(mySeurat, features=c("Slc17a6"), idents=c("21","26"), pt.size=0, cols=c("gray","gray"))

glutfdat <- glutf$data
pdf(file=paste0(output,"femaleglutcust.pdf"))
ggplot(glutfdat, aes(x=Slc17a6, y=factor(ident, levels=rev(levels(factor(ident)))))) + geom_violin(trim=TRUE, fill="gray", linetype="blank") + theme_classic() +theme(axis.text.x=element_text(size=20), axis.title.y=element_blank()) + xlim(0,2.5)
dev.off()
write.csv(glutfdat, file=paste0(output, "glutfvlndata.csv"))
pdf(file=paste0(output,"Kiss1black.pdf"))
VlnPlot(mySeurat, features=c("Kiss1"), idents=c("26"), pt.size=0, cols=c("black"))
dev.off()
pdf(file=paste0(output,"Ror1black.pdf"))
VlnPlot(mySeurat, features=c("Ror1"), idents=c("10"), pt.size=0, cols=c("black"))
dev.off()
pdf(file=paste0(output,"femaleupSlc18a2.pdf"))
VlnPlot(mySeurat, features=c("Slc18a2"), idents=c("21","26"), pt.size=0, cols=c("gray","gray"))
dev.off()
pdf(file=paste0(output,"labeledUMAPFigS4.pdf"))
DimPlot(mySeurat, label=TRUE, reduction="umap",label.size=6)
dev.off()
pdf(file=paste0(output,"inhibexcitcluster.pdf"))
DimPlot(mySeurat, cols=c("0"="blue", "1"="red", "2"="blue", "3"="red", "4"="red", "5"="red", "6"="red", "7"="blue", "8"="blue", "9"="blue", "10"="red", "11"="red", "12"="red", "13"="red", "14"="red", "15"="blue", "16"="blue", "17"="red", "18"="blue", "19"="red", "20"="red", "21"="red", "22"="blue", "23"="red", "24"="red", "25"="red", "26"="red", "27"="red", "28"="blue", "29"="red", "30"="blue", "31"="red", "32"="red", "33"="red", "34"="red", "35"="red","36"="blue", "37"="blue", "38"="blue"), label=FALSE, reduction="umap")
dev.off() 
genelist <- read.table("Genelists/POA_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
MvP <- read.table("Genelists/POA_MvP_1.5cutoff.txt", header=FALSE)
MvU <- read.table("Genelists/POA_MvU_1.5cutoff.txt", header=FALSE)
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
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/AllDEGsavg.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/POA/AllDEGsavg.txt", header=TRUE, row.names=1)
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
