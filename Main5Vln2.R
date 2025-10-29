#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_CelltypesDraftArcfiltered_regressed_reorderedJune"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)
library(patchwork)


BNSTgenes <- c("Rbm20","Plekhg1","Prok2","Cck","Phf21b","Slco2a1","Tac1","Asic4","Kank1","Sst","Abca1","Amot","Fam107b","Prkcd","Sertm1","Tns1","Sytl4","Npr3","Plscr4")
MeAgenes <- c("Sst","Map3k15","Tspan18","Vgll3","Efemp1")
POAgenes <- c("Efemp1","Mob3b","Moxd1","St3gal1","Thbs4","Mrc2","Ttn","Neb","Kiss1","Kl","Slc17a8","Abtb2","Penk")
VMHgenes <- c("Cckar","Tnfaip8l3","Mme","Cgnl1","Popdc3","Trim36")


BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(BNST)
BNST <- RenameIdents(BNST, new.cluster.ids)
DefaultAssay(BNST) <- "RNA"
BNST <- NormalizeData(BNST)

counts <- BNST@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,BNSTgenes]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
BNSTout <- cbind(BNSTgenes,max.df)
BNSTout$CPM <- "CPM"

BNSTcpm <- ggplot(BNSTout, aes(y=factor(BNSTgenes, levels=BNSTgenes), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=15) + scale_color_viridis(limits=c(2.3,3.6)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

pdf(file="Seurat/BNSTSpecDEGcpm.pdf")
ggplot(BNSTout, aes(y=factor(BNSTgenes, levels=BNSTgenes), x=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=20) + scale_color_viridis() +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()


plot <- VlnPlot(BNST, features=rev(BNSTgenes), pt.size=0, ncol=1)
plotf <- plot & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & geom_violin(trim=TRUE, fill="red", scale="width") 

pdf(file="Seurat/BNST_SpecDEGVLnnew.pdf", width=37, height=30)
(plotf|BNSTcpm)  + plot_layout(widths=c("1","0.05")) 
dev.off()


VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(VMH)
VMH <- RenameIdents(VMH, new.cluster.ids)
DefaultAssay(VMH) <- "RNA"
VMH <- NormalizeData(VMH)

counts <- VMH@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,VMHgenes]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
vmhout <- cbind(VMHgenes,max.df)

vmhout$CPM <- "CPM"

VMHcpm <- ggplot(vmhout, aes(y=factor(VMHgenes, levels=VMHgenes), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=15) + scale_color_viridis(limits=c(2.3,3.6)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

pdf(file="Seurat/VMHSpecDEGcpm.pdf")
ggplot(vmhout, aes(y=factor(VMHgenes, levels=VMHgenes), x=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=20) + scale_color_viridis() +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()

plot <- VlnPlot(VMH, features=rev(VMHgenes), pt.size=0, ncol=1)
plotf <- plot & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & geom_violin(trim=TRUE, fill="red", scale="width", weight=0.5)

pdf(file="Seurat/VMH_SpecDEGVLnnew.pdf", width=37, height=9.5)
(plotf|VMHcpm)  + plot_layout(widths=c("1","0.05")) 
dev.off()

MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
names(new.cluster.ids) <- levels(MeA)
MeA <- RenameIdents(MeA, new.cluster.ids)
DefaultAssay(MeA) <- "RNA"
MeA <- NormalizeData(MeA)

counts <- MeA@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,MeAgenes]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
meaout <- cbind(MeAgenes,max.df)

meaout$CPM <- "CPM"

MeAcpm <- ggplot(meaout, aes(y=factor(MeAgenes, levels=MeAgenes), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=15) + scale_color_viridis(limits=c(2.3,3.6)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

pdf(file="Seurat/MeASpecDECPM.pdf")
ggplot(meaout, aes(y=factor(MeAgenes, levels=MeAgenes), x=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=20) + scale_color_viridis() +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()


plot <- VlnPlot(MeA, features=rev(MeAgenes), pt.size=0, ncol=1)
plotf <- plot & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & geom_violin(trim=TRUE, fill="red", scale="width", weight=0.25)

pdf(file="Seurat/MeA_SpecDEGVLnnew.pdf", width=37, height=7.8)
(plotf|MeAcpm)  + plot_layout(widths=c("1","0.05")) 
dev.off()

POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(POA)
POA <- RenameIdents(POA, new.cluster.ids)
DefaultAssay(POA) <- "RNA"
POA <- NormalizeData(POA)

counts <- POA@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,POAgenes]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
POAout <- cbind(POAgenes,max.df)

POAout$CPM <- "CPM"

POAcpm <- ggplot(POAout, aes(y=factor(POAgenes, levels=POAgenes), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=15) + scale_color_viridis(limits=c(2.3,3.6)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

pdf(file="Seurat/POA_SpecDEGcpm.pdf")
ggplot(POAout, aes(y=factor(POAgenes, levels=POAgenes), x=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=20) + scale_color_viridis() +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()

plot <- VlnPlot(POA, features=rev(POAgenes), pt.size=0, ncol=1)
plotf <- plot & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & geom_violin(trim=TRUE, fill="red", scale="width", weight=0.25)

pdf(file="Seurat/POA_SpecDEGVLnnew.pdf", width=37, height=20.5)
(plotf|POAcpm)  + plot_layout(widths=c("1","0.05")) 
dev.off()