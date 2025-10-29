#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_FigS7Main7_final"

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

DefaultAssay(mySeurat) <- "RNA"

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
mySeurat <- NormalizeData(mySeurat)

edegmarkers <- read.table("Genelists/Cckar_eDEGmarkers.txt", header=FALSE)
edegmarkers <- unlist(edegmarkers)
pdf(file=paste0(output,"Tonsofvagplots.pdf"), width=20,height=300)
VlnPlot(mySeurat, features=edegmarkers,ncol=1,pt.size=0)
dev.off()

Male <- subset(mySeurat, subset=sex=="Male")
Female <- subset(mySeurat, subset=sex=="Female")

Mc <- FeaturePlot(Male, features=c("Cckar"))
mcdat <- Mc$data
Fc <- FeaturePlot(Female, features=c("Cckar"))
fcdat <- Fc$data

p1 <- ggplot(mcdat, aes(x=UMAP_1, y=UMAP_2, color=Cckar)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("light gray","red"), limits=c(0, 1)) + cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank())


p2 <- ggplot(fcdat, aes(x=UMAP_1, y=UMAP_2, color=Cckar)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("light gray","red"), limits=c(0, 1)) + cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank())

pdf(file=paste0(output,"UmapCckar.pdf"), width=12)
p1 + p2
dev.off()



v <- VlnPlot(mySeurat, features=c("Cckar"), pt.size=0)
vat <- v$data


vv1 <- ggplot(vat, aes(x=Cckar, y=factor(ident, levels=mylevels), fill=ident) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("deeppink1",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1)) +theme(legend.position="none")

pdf(file=paste0(output,"vertvlntest.pdf"), height=15)
vv1
dev.off()
vlnlevels <- c(6,20,25)

v <- VlnPlot(mySeurat, features=c("Abtb2"), pt.size=0)
vat <- v$data



vv2 <- ggplot(vat, aes(x=Abtb2, y=factor(ident, levels=mylevels), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("deeppink1",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

v <- VlnPlot(mySeurat, features=c("Myo1h"), pt.size=0)
vat <- v$data


vv3 <- ggplot(vat, aes(x=Myo1h, y=factor(ident, levels=mylevels), fill=ident)) + geom_violin(trim=TRUE, scale="width") +scale_fill_manual(values=c(rep("dark gray", 5), rep("deeppink1",1), rep("dark gray",13), rep("deeppink1",1), rep("dark gray",4), rep("green", 1), rep("dark gray",2)))  + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_x_continuous(breaks=seq(0,4,1))+theme(legend.position="none")

pdf(file=paste0(output,"vertvln.pdf"), height=15)
vv1
dev.off()


pdf(file=paste0(output,"vertvlntest.pdf"), height=15)
vv1|vv2|vv3
dev.off()


v <- VlnPlot(mySeurat, features=c("Cckar"), pt.size=0)
vat <- v$data



fv1 <- ggplot(vat, aes(y=Cckar, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


v <- VlnPlot(mySeurat, features=c("Arg2"), pt.size=0)
vat <- v$data


fv2 <- ggplot(vat, aes(y=Arg2, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Tmem215"), pt.size=0)
vat <- v$data


fv3 <- ggplot(vat, aes(y=Tmem215, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Gfra1"), pt.size=0)
vat <- v$data


fv4 <- ggplot(vat, aes(y=Gfra1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Mylk"), pt.size=0)
vat <- v$data


fv5 <- ggplot(vat, aes(y=Mylk, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Phf21b"), pt.size=0)
vat <- v$data


fv6 <- ggplot(vat, aes(y=Phf21b, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
v <- VlnPlot(mySeurat, features=c("Nfil3"), pt.size=0)
vat <- v$data


fv7 <- ggplot(vat, aes(y=Nfil3, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Tmod1"), pt.size=0)
vat <- v$data


fv8 <- ggplot(vat, aes(y=Tmod1, x=factor(ident, levels=mylevels))) + geom_violin(trim=TRUE, fill="dark gray", scale="width") + theme_classic() +theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"S1_eDEGsstacked.pdf"), width=30, height=20)
fv1/fv2/fv3/fv4/fv5/fv6/fv7/fv8
dev.off()


v <- VlnPlot(mySeurat, features=c("Cckar"), idents=c("6","20","25"))

vat <- v$data

vm1 <- ggplot(vat, aes(y=Cckar, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


v <- VlnPlot(mySeurat, features=c("Abtb2"), idents=c("6","20","25"))

vat <- v$data

vm2 <- ggplot(vat, aes(y=Abtb2, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Cyp19a1"), idents=c("6","20","25"))

vat <- v$data

vm3 <- ggplot(vat, aes(y=Cyp19a1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Myo1h"), idents=c("6","20","25"))

vat <- v$data

vm4 <- ggplot(vat, aes(y=Myo1h, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
blank <- ggplot() + theme_void()
pdf(file=paste0(output,"MarkerVlnStack.pdf"), height=13)
vm1/vm2/vm4/blank
dev.off()


v <- VlnPlot(mySeurat, features=c("Arg2"), idents=c("6","20","25"))

vat <- v$data

ve1 <- ggplot(vat, aes(y=Arg2, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))


v <- VlnPlot(mySeurat, features=c("Tmem215"), idents=c("6","20","25"))

vat <- v$data

ve2 <- ggplot(vat, aes(y=Tmem215, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Gfra1"), idents=c("6","20","25"))

vat <- v$data

ve3 <- ggplot(vat, aes(y=Gfra1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Mylk"), idents=c("6","20","25"))

vat <- v$data

ve4 <- ggplot(vat, aes(y=Mylk, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
pdf(file=paste0(output,"FrUP1.pdf"), height=13)
ve1/ve2/ve3/ve4
dev.off()

v <- VlnPlot(mySeurat, features=c("Phf21b"), idents=c("6","20","25"))

vat <- v$data

ve5 <- ggplot(vat, aes(y=Phf21b, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Nfil3"), idents=c("6","20","25"))

vat <- v$data

ve6 <- ggplot(vat, aes(y=Nfil3, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))
v <- VlnPlot(mySeurat, features=c("Tmod1"), idents=c("6","20","25"))

vat <- v$data

ve7 <- ggplot(vat, aes(y=Tmod1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Tmod1"), idents=c("6","20","25"))

vat <- v$data

ve8 <- ggplot(vat, aes(y=Tmod1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"FrUP2.pdf"), height=13)
ve5/ve6/ve7/blank
dev.off()

v <- VlnPlot(mySeurat, features=c("Cgnl1"), idents=c("6","20","25"))

vat <- v$data

vd1 <- ggplot(vat, aes(y=Cgnl1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Popdc3"), idents=c("6","20","25"))

vat <- v$data

vd2 <- ggplot(vat, aes(y=Popdc3, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Trim36"), idents=c("6","20","25"))

vat <- v$data

vd3 <- ggplot(vat, aes(y=Trim36, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

v <- VlnPlot(mySeurat, features=c("Klhl1"), idents=c("6","20","25"))

vat <- v$data

vd4 <- ggplot(vat, aes(y=Klhl1, x=factor(ident, levels=vlnlevels), fill=ident)) + geom_violin(trim=TRUE,  scale="width") +scale_fill_manual(values=c("deeppink1","deeppink1","green"))+ theme_classic() +theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())+theme(legend.position="none") +theme(axis.line=element_line(color="black", size=1)) + scale_y_continuous(breaks=seq(0,4,1))

pdf(file=paste0(output,"FuUP1.pdf"), width=13)
vd1/vd2/vd3/vd4
dev.off()