#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_Fig7main_recluster"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)
library(tidyr)
library(tidyverse)
library(patchwork)
mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
dat <- VlnPlot(mySeurat, features=c("nCount_RNA"))
data <- dat$data
data$val <- "0"


p1 <- ggplot(data, aes(y=nCount_RNA, x=val)) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_text(size=12), axis.text.y=element_text(size=60), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=2))

dat2 <- VlnPlot(mySeurat, features=c("nFeature_RNA"))
data2 <- dat2$data
data2$val <- "0"

p2 <- ggplot(data2, aes(y=nFeature_RNA, x=val)) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
pdf(file=paste0(output,"genevln.pdf"))
p1 + p2
dev.off()

mySeurat <- NormalizeData(mySeurat)

pdf(file=paste0(output,"oxtrvln.pdf"), width=30)
VlnPlot(mySeurat, features=c("Oxtr"), pt.size=1)
dev.off()
pdf(file=paste0(output,"oxtrvlnsplit.pdf"), width=30)
VlnPlot(mySeurat, features=c("Oxtr"), pt.size=1, split.by="Hormone")
dev.off()
pdf(file=paste0(output,"oxtrftr.pdf"))
FeaturePlot(mySeurat, features=c("Oxtr"), order=TRUE, cols=c("light blue", "red"))
dev.off()
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
pdf(file=paste0(output,"sexbias.pdf"))
DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="light gray", "4"="light gray", "5"="light gray", "6"="light gray", 
	"7"="light gray", "8"="deeppink1", "9"="light gray", "10"="light gray", "11"="light gray", "12"="light gray", "13"="light gray", 
	"14"="light gray", "15"="light gray", "16"="light gray", "17"="light gray", "18"="light gray", "19"="blue", "20"="light gray", 
	"21"="deeppink1", "22"="light gray", "23"="light gray", "24"="light gray", "25"="light gray", "26"="green", "27"="light gray"), label=FALSE, reduction="umap")
dev.off() 
pdf(file=paste0(output,"sexbiasaltMale.pdf"))
DimPlot(mySeurat, group.by="Hormone", shuffle=TRUE, cols=c("Blue","Gray","Gray"), label=FALSE)
dev.off()
pdf(file=paste0(output,"sexbiasaltPrimed.pdf"))
DimPlot(mySeurat, group.by="Hormone", shuffle=TRUE, cols=c("gray","deeppink1","Gray"), label=FALSE)
dev.off()
pdf(file=paste0(output,"sexbiasaltUnprimed.pdf"))
DimPlot(mySeurat, group.by="Hormone", shuffle=TRUE, cols=c("gray","gray","Green"), label=FALSE)
dev.off()
MaleUP <- read.table("DEGmodules/VMH_MaleupvPrimed.txt", header=FALSE)
MaleUP <- unlist(MaleUP)
MaleUP
maleUP.filtered <- intersect(MaleUP,genes.10x)
maleUP.filtered
FrUP <- read.table("DEGmodules/VMH_UpPrimedvFu.txt", header=FALSE)
FrUP <- unlist(FrUP)
FrUP.filtered <- intersect(FrUP, genes.10x)
FuUP <- read.table("DEGmodules/VMH_UpFuvPrimed.txt", header=FALSE)
FuUP <- unlist(FuUP)
FuUP.filtered <- intersect(FuUP, genes.10x)

mySeurat <- AddModuleScore(mySeurat, features=list(maleUP.filtered), assay="RNA", name="maleUP")
mySeurat <- AddModuleScore(mySeurat, features=list(FrUP.filtered), assay="RNA", name="FrUP")
mySeurat <- AddModuleScore(mySeurat, features=list(FuUP.filtered), assay="RNA", name="FuUP")
pdf(file=paste0(output,"modvln.pdf"), width=30, height=20)
VlnPlot(mySeurat, features=c("maleUP1","FrUP1","FuUP1"), pt.size=0, ncol=1)
dev.off()
pdf(file=paste0(output,"modftr.pdf"), width=12, height=20)
FeaturePlot(mySeurat, features=c("maleUP1","FrUP1","FuUP1"), cols=c("light blue", "red"))
dev.off()

mySeurat[['modulescores']] <- CreateAssayObject(data = t(x = FetchData(object = mySeurat, vars = c("maleUP1","FrUP1","FuUP1"))))
pdf(file=paste0(output,"Moduleheatmaptest.pdf"))
DoHeatmap(object = mySeurat, features = c("maleUP1","FrUP1","FuUP1"), assay = 'modulescores', slot = 'data')
dev.off()


data <- mySeurat@assays[["modulescores"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)
mods <- c("maleUP1","FrUP1","FuUP1")

avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=mods, .funs=c("mean"))
avg
avg <- avg %>% remove_rownames %>% column_to_rownames (var="clusters")
avg
avg <- scale(avg)
write.csv(avg, file="Heatmaps/Seurat/VMH/ScaledModscore.csv")
pdf(file=paste0(output,"scaledmodheatmap.pdf"))
pheatmap(avg, color=viridis(2000), cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mylevels <- c(1,2,3,5,6,10,11,14,15,16,18,19,20,23,24,25,26,0,4,7,8,9,12,13,17,21,22)

glut <- VlnPlot(mySeurat, features=c("Esr1"), pt.size=0)
glutdat <- glut$data


mv1 <- ggplot(glutdat, aes(y=Esr1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_text(size=80, color="black"), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,4,1))


glut <- VlnPlot(mySeurat, features=c("Pgr"), pt.size=0)
glutdat <- glut$data


mv2 <- ggplot(glutdat, aes(y=Pgr, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_text(size=80, color="black"), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,3,1))


glut <- VlnPlot(mySeurat, features=c("Gad1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Gad1vlndata.csv"))
pdf(file=paste0(output,"Gad1vln.pdf"), width=40)
ggplot(glutdat, aes(y=Gad1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_text(size=30), axis.title.y=element_blank())+theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Cyp19a1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Cyp19a1vlndata.csv"))
pdf(file=paste0(output,"Cyp19a1vln.pdf"), width=40)
ggplot(glutdat, aes(y=Cyp19a1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30),axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Slc17a6"), pt.size=0)
glutdat <- glut$data

mv3 <- ggplot(glutdat, aes(y=Slc17a6, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=80, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,4,1))

glut <- VlnPlot(mySeurat, features=c("Phf21b"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Phf21bfullvlndata.csv"))

v1 <- ggplot(glutdat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2.5,1))


glut <- VlnPlot(mySeurat, features=c("Samd9l"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Samd9lfullvlndata.csv"))

v2 <- ggplot(glutdat, aes(y=Samd9l, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2.2,1))


glut <- VlnPlot(mySeurat, features=c("Mylk"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Mylkfullvlndata.csv"))

v3 <- ggplot(glutdat, aes(y=Mylk, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2.5,1))


glut <- VlnPlot(mySeurat, features=c("Tmem215"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Tmem215fullvlndata.csv"))

v4 <- ggplot(glutdat, aes(y=Tmem215, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,1.5,1))

glut <- VlnPlot(mySeurat, features=c("St3gal1"), pt.size=0)
glutdat <- glut$data
v5 <- ggplot(glutdat, aes(y=St3gal1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2,1))

glut <- VlnPlot(mySeurat, features=c("Ptp4a1"), pt.size=0)
glutdat <- glut$data
v6 <- ggplot(glutdat, aes(y=Ptp4a1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,2,1))

glut <- VlnPlot(mySeurat, features=c("Mical2"), pt.size=0)
glutdat <- glut$data
v7 <- ggplot(glutdat, aes(y=Mical2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,3,1))

glut <- VlnPlot(mySeurat, features=c("Gfra1"), pt.size=0)
glutdat <- glut$data
v8 <- ggplot(glutdat, aes(y=Gfra1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=45, color="black"),axis.title.y=element_blank())  +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,3,1))


pdf(file=paste0(output,"eDEGvlnsstacked.pdf"), width=60, height=35)
(v1/v2/v3/v4/v5/v6/v7/v8)

dev.off()
glut <- VlnPlot(mySeurat, features=c("Cckar"), cols=c("dark gray"), pt.size=0)

glutdat <- glut$data

mv4 <- ggplot(glutdat, aes(y=Cckar, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=80, color="black"), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=2)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"Marker_vln_stacked.pdf"), width=80, height=30)
(mv1/mv3/mv2/mv4)
dev.off()

glut <- VlnPlot(mySeurat, features=c("Rprm"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Rprmvlndata.csv"))
pdf(file=paste0(output,"Rprmvln.pdf"), width=40)
ggplot(glutdat, aes(y=Rprm, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()


glut <- VlnPlot(mySeurat, features=c("Iigp1"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Iigp1vlndata.csv"))
pdf(file=paste0(output,"Iigp1vln.pdf"), width=40)
ggplot(glutdat, aes(y=Iigp1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Phf21b"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Phf21bvlndata.csv"))
pdf(file=paste0(output,"Phf21bvln.pdf"), width=40)
ggplot(glutdat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray",  scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Zfp804a"), pt.size=0)
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Zfp804avlndata.csv"))
pdf(file=paste0(output,"Zfp804avln.pdf"), width=40)
ggplot(glutdat, aes(y=Zfp804a, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()



glut <- VlnPlot(mySeurat, features=c("Cckar"), idents=c("8","21","26"), pt.size=0, ncol=1, cols=c("dark gray","dark gray","dark gray"))
pdf(file=paste0(output,"CZIcckarvlnraw.pdf"))
glut
dev.off()
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "cckarCZIvlndata.csv"))
pdf(file=paste0(output,"CckarCZIvln.pdf"))
ggplot(glutdat, aes(y=Cckar, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()


glut <- VlnPlot(mySeurat, features=c("Myo1h"), idents=c("8","21","26"), pt.size=0, ncol=1, cols=c("dark gray","dark gray","dark gray"))
pdf(file=paste0(output,"CZImyo1hvlnraw.pdf"))
glut
dev.off()
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "myo1hCZIvlndata.csv"))
pdf(file=paste0(output,"Myo1hCZIvln.pdf"))
ggplot(glutdat, aes(y=Myo1h, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

glut <- VlnPlot(mySeurat, features=c("Phf21b"), idents=c("8","21","26"), pt.size=0, ncol=1, cols=c("dark gray","dark gray","dark gray"))
pdf(file=paste0(output,"Phf21bCZIvlnraw.pdf"))
glut
dev.off()
glutdat <- glut$data
write.csv(glutdat, file=paste0(output, "Phf21bCZIvlndata.csv"))
pdf(file=paste0(output,"Phf21bCZIvln.pdf"))
ggplot(glutdat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_text(size=12), axis.text.y=element_text(size=30), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1))
dev.off()

AvgAllDEGs <- AverageExpression(mySeurat, features=c("Cckar","Rprm","Rnf19a","Iigp1","Maml2","Mc4r","Mad2l1","Tmod1"), assays=c("RNA"))
AvgAllDEGs
write.csv(AvgAllDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/Cckarmarkerssavg.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/Cckarmarkerssavg.txt", header=TRUE, row.names=1)
data <- as.matrix(data)
dim(data)
head(data)
data <- scale(t(data))
write.csv(data, file=paste0(output,"ScaledCckarMarkers.csv"))
idents <- mySeurat@active.ident
idents <- as.matrix(idents)
head(idents)
pdf(file=paste0(output,"Cckarheatmaptrowclustonly.pdf"))
pheatmap(data, cluster_cols=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"Cckarheatmapnoclust.pdf"))
pheatmap(data, cluster_cols=FALSE, cluster_rows=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"Cckarheatmapt.pdf"), width=30, height=30)
pheatmap(data, border_color=NA, color=viridis(500))
dev.off()






AvgselectDEGs <- AverageExpression(mySeurat, features=c("Mad2l1","Stat5a","Foxo1","Tmem215","Tnxb","Ezr","Cyp19a1","Gldn","Moxd1","Gal","Pdyn","Il1rap","Rbms3","Cartpt","Calb1","Greb1","Zfp804a","Klhl1","Nxph1","Cdh10","B3galt1","Cab39l","Cgnl1","Popdc3","Trim36","Iqschfp"), assays=c("RNA"))
AvgselectDEGs
write.csv(AvgselectDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/SelectDEGsavg2.txt")
data <- read.csv(file="/scratch/users/knoedler/Heatmaps/Seurat/VMH/SelectDEGsavg2.txt", header=TRUE, row.names=1)
data <- as.matrix(data)
dim(data)
head(data)
data <- scale(t(data))
pdf(file=paste0(output,"SelectDEGsavgheatmaptrowclustonly.pdf"))
pheatmap(data, cluster_cols=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"SelectDEGsavgheatmapnoclust.pdf"))
pheatmap(data, cluster_cols=FALSE, cluster_rows=FALSE, border_color=NA, color=viridis(500))
dev.off()
pdf(file=paste0(output,"SelectDEGsavgheatmapt.pdf"), width=30, height=30)
pheatmap(data, border_color=NA, color=viridis(500))
dev.off()