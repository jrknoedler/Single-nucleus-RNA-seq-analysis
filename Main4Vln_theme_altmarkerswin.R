#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "X:/Seurat/VMH_IndependentAnalysis/VMH_CelltypesDraftArcfiltered_regressed_reorderedJune"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(matrixStats)
memory.limit(64000)
blank <- ggplot() + theme_void() +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) 
phenotypes <- c("Slc17a6","Slc17a7","Gad1","Cyp19a1","Esr2","Th")
phenotypes <- unlist(phenotypes)
phenotypes <- rev(phenotypes)
BNST <- readRDS("X:/Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(BNST)
BNST <- RenameIdents(BNST, new.cluster.ids)
DefaultAssay(BNST) <- "RNA"
BNST <- NormalizeData(BNST)


#BNSTtypemarkers <- read.table("X:/DotPlotMarkerLists/BNST_ReorderedFinal.txt", header=FALSE)
BNSTtypemarkers <- read.table("X:/DotPlotMarkerLists/BNST_ReorderedFinalalternate_CPM.txt", header=FALSE)
BNSTtypemarkers <- unlist(BNSTtypemarkers)
BNSTtypemarkers.m <- as.matrix(BNSTtypemarkers)

BNSTtypemarkersnt <- read.table("X:/DotPlotMarkerLists/BNST_ReorderedFinal_CPM_excitmarkers.txt", header=FALSE)
BNSTtypemarkersnt <- unlist(BNSTtypemarkersnt)



BNSTlevels <- c(6,15,20,36,1,2,3,4,5,7,8,9,10,11,12,13,14,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)
#BNSTlevels <- rev(BNSTlevels)

#BNST@ident <- factor(x=BNST@ident, levels=BNSTlevels)
Idents(BNST) <- factor(Idents(BNST), levels=BNSTlevels)


counts <- BNST@assays[["RNA"]]@counts
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- as.matrix(CPM)
CPM <- t(CPM)

genes.m <- CPM[,BNSTtypemarkers]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
BNSTout <- cbind(BNSTtypemarkers,max.df)

genes.ph <- CPM[,phenotypes]
genes.ph <- as.matrix(genes.ph)
max <- colMaxs(genes.ph)
max.df <- data.frame(max)
BNSTph <- cbind(phenotypes,max.df)
BNSTph$CPM <- "CPM"

BNSTphCPM <- ggplot(BNSTph, aes(y=factor(phenotypes, levels=rev(phenotypes)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(2,3.5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())


plotph <- VlnPlot(BNST, features=phenotypes, pt.size=0, ncol=1)
plotphf <- plotph & theme_nothing() & theme(axis.line.x=element_line(color="black", size=1)) & scale_fill_manual(values=c(rep("royalblue3",4), rep("red",32)))

pdf(file="X:/Seurat/BNST_Phenotypes.pdf", width=50, height=20)
(plotphf|BNSTphCPM)  + plot_layout(widths=c("1","0.05")) 
dev.off()


BNSTout$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
BNSTcpm <- ggplot(BNSTout, aes(y=factor(BNSTtypemarkers, levels=BNSTtypemarkers), x=CPM, color=log10(max))) +theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(1,4.5)) +  theme(axis.line  = element_blank())  +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) + theme(axis.title=element_blank()) +  theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

genes.mnt <- CPM[,BNSTtypemarkersnt]

genes.fnt <- as.matrix(genes.mnt)
maxnt <- colMaxs(genes.fnt)
max.dfnt <- data.frame(maxnt)
BNSToutnt <- cbind(BNSTtypemarkersnt,max.dfnt)

BNSToutnt$CPM <- "CPM"

plot <- VlnPlot(BNST, features=rev(BNSTtypemarkers), pt.size=0, ncol=1)
plotf <- plot & theme_nothing() & theme(axis.line.x=element_line(color="black", size=1)) & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) & scale_fill_manual(values=c(rep("royalblue3",4), rep("red",32)))
#BNSTvln <- (plotf|BNSTcpm) & theme(panel.background=element_rect(fill="transparent"))+ plot_layout(widths=c("1","0.05")) 
bnstc <- (plotf/blank) + plot_layout(heights=c("1","0.02"))
#ggsave("X:/Seurat/BNST_ViolinTest.pdf",BNSTvln, bg='transparent')

pdf(file="X:/Seurat/BNST_ViolinTest_alternate.pdf", width=60, height=163)
(plotf|BNSTcpm) + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()

pdf(file="X:/Seurat/BNST_ViolinTest_alternate.pdf", width=57, height=163)
plotf & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()


pdf(file="X:/Seurat/BNST_Violin_altCPM", width=3, height=167)
BNSTcpm
dev.off()

pdf(file="X:/Seurat/BNST_ViolinTest_alternatenames.pdf", width=60, height=163)
(plotf|BNSTcpm) + plot_layout(widths=c("1","0.05")) & cowplot::theme_cowplot()
dev.off()

BNSTcpm <- ggplot(BNSTout, aes(y=factor(BNSTtypemarkers, levels=BNSTtypemarkers), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(1,4.5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
BNSTcpmnt <- ggplot(BNSToutnt, aes(y=factor(BNSTtypemarkersnt, levels=BNSTtypemarkersnt), x=CPM, color=log10(maxnt))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(2.3,4)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

plot2 <- VlnPlot(BNST, features=BNSTtypemarkersnt, pt.size=0, ncol=1)
plotf2 <- plot2 & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & scale_fill_manual(values=c(rep("royalblue3",4), rep("red",32)))

pdf(file="X:/Seurat/BNST_ViolinEdittestwide.pdf", width=95, height=38.5)
(plotf|BNSTcpm)  + plot_layout(widths=c("1","0.05")) 
dev.off()

pdf(file="X:/Seurat/BNSTCPMwidescale.pdf", height=78)
ggplot(BNSToutnt, aes(y=factor(BNSTtypemarkersnt, levels=BNSTtypemarkersnt), x=CPM, color=log10(maxnt))) +cowplot::theme_cowplot() + geom_point(size=30) + scale_color_viridis(limits=c(2.3,4)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()
VMH <- readRDS("X:/Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(VMH)
VMH <- RenameIdents(VMH, new.cluster.ids)
DefaultAssay(VMH) <- "RNA"
VMH <- NormalizeData(VMH)
#vmhtypemarkers <- read.table("X:/DotPlotMarkerLists/VMH_reorderedfinal.txt", header=FALSE)
vmhtypemarkers <- read.table("X:/DotPlotMarkerLists/VMH_reorderedfinaledited_CPM.txt", header=FALSE)
vmhtypemarkers <- unlist(vmhtypemarkers)
vmhtypemarkers <- unlist(vmhtypemarkers)
vmhtypemarkers.m <- as.matrix(vmhtypemarkers)

vmhtypemarkersnt <- read.table("X:/DotPlotMarkerLists/VMH_reorderedfinal_CPM_excitmarkers.txt", header=FALSE)
vmhtypemarkersnt <- unlist(vmhtypemarkersnt)



VMHlevels <- c(2,3,4,6,7,11,12,15,16,17,19,20,21,24,25,26,27,1,5,8,9,10,13,14,18,22,23)
#VMHlevels <- rev(VMHlevels)
#VMH@ident <- factor(x=VMH@ident, levels=VMHlevels)
Idents(VMH) <- factor(Idents(VMH), levels=VMHlevels)
pdf(file="X:/Seurat/VMH_ViolinTest_altmarkers.pdf", width=60, height=126)
VlnPlot(VMH, features=vmhtypemarkers, pt.size=0, ncol=1)
dev.off()

plot <- VlnPlot(VMH, features=rev(vmhtypemarkers), pt.size=0, ncol=1)
plotf <- plot& theme_nothing() & theme(axis.line.x=element_line(color="black", size=1)) & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) &  scale_fill_manual(values=c(rep("royalblue3",17), rep("red",10)))


counts <- VMH@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- as.matrix(CPM)
CPM <- t(CPM)

genes.m <- CPM[,vmhtypemarkers]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
vmhout <- cbind(vmhtypemarkers,max.df)

vmhout$CPM <- "CPM"

genes.mnt <- CPM[,vmhtypemarkersnt]

genes.fnt <- as.matrix(genes.mnt)
maxnt <- colMaxs(genes.fnt)
max.dfnt <- data.frame(maxnt)
vmhoutnt <- cbind(vmhtypemarkersnt,max.dfnt)

vmhoutnt$CPM <- "CPM"


#typemarkers <- typemarkers[,data$features.plot[["order"]]]
vmhcpm <- ggplot(vmhout, aes(y=factor(vmhtypemarkers, levels=vmhtypemarkers), x=CPM, color=log10(max))) +theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(1,4.5)) +  theme(axis.line  = element_blank())  +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) + theme(axis.title=element_blank()) +  theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())


pdf(file="X:/Seurat/VMH_ViolinEdittest_alternate.pdf", width=60, height=125.6)
(plotf|vmhcpm) + plot_layout(widths=c("1","0.05")) & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()
#pdf(file="X:/Seurat/VMH_ViolinEdittest_alternate.pdf", width=57, height=125.6)
#plotf
#dev.off()
pdf(file="X:/Seurat/VMH_CPMalt.pdf",width=5, height=130)
vmhcpm & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()


#pdf(file="X:/Seurat/VMH_ViolinEdittest_names.pdf", width=60, height=125.6)
#(plotf|vmhcpm) + plot_layout(widths=c("1","0.05")) & cowplot::theme_cowplot()
#dev.off()

genes.ph <- CPM[,phenotypes]
genes.ph <- as.matrix(genes.ph)
max <- colMaxs(genes.ph)
max.df <- data.frame(max)
VMHph <- cbind(phenotypes,max.df)
VMHph$CPM <- "CPM"

VMHphCPM <- ggplot(VMHph, aes(y=factor(phenotypes, levels=rev(phenotypes)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(2,3.5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())


plotph <- VlnPlot(VMH, features=phenotypes, pt.size=0, ncol=1)
plotphf <- plotph & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & scale_fill_manual(values=c(rep("royalblue3",17), rep("red",10)))

pdf(file="X:/Seurat/VMH_Phenotypes.pdf", width=50, height=20)
(plotphf|VMHphCPM) + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()


plot2 <- VlnPlot(VMH, features=vmhtypemarkersnt, pt.size=0, ncol=1)
plotf2 <- plot2 & theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & scale_fill_manual(values=c(rep("royalblue3",17), rep("red",10)))
vmhcpm <- ggplot(vmhout, aes(y=factor(vmhtypemarkers, levels=vmhtypemarkers), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(1,4.5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
vmhcpmnt <- ggplot(vmhoutnt, aes(y=factor(vmhtypemarkersnt, levels=vmhtypemarkersnt), x=CPM, color=log10(maxnt))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(2.3,4)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

pdf(file="X:/Seurat/VMH_ViolinEdittestwid.pdf", width=95, height=30)
(plotf|vmhcpm)  + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()
pdf(file="X:/Seurat/VMHCPMwidescale.pdf", height=78)
ggplot(vmhoutnt, aes(y=factor(vmhtypemarkersnt, levels=vmhtypemarkersnt), x=CPM, color=log10(maxnt))) +cowplot::theme_cowplot() + geom_point(size=30) + scale_color_viridis() +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())
dev.off()

MeA <- readRDS("X:/Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
names(new.cluster.ids) <- levels(MeA)
MeA <- RenameIdents(MeA, new.cluster.ids)
DefaultAssay(MeA) <- "RNA"
MeA <- NormalizeData(MeA)
typemarkers <- read.table("X:/DotPlotMarkerLists/MeA_ReorderedFinalalternate_CPM.txt", header=FALSE)
typemarkers <- unlist(typemarkers)
typemarkers <- unlist(typemarkers)
typemarkers.m <- as.matrix(typemarkers)
MeAlevels <- c(1,4,7,9,10,13,14,15,16,18,21,22,23,24,29,30,31,32,33,2,3,5,6,8,11,12,17,19,20,25,26,27,28,34)
#MeAlevels <- rev(MeAlevels)
Idents(MeA) <- factor(Idents(MeA), levels=MeAlevels)
plot <- VlnPlot(MeA, features=rev(typemarkers), pt.size=0, ncol=1)


counts <- MeA@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- as.matrix(CPM)
CPM <- t(CPM)

genes.m <- CPM[,typemarkers]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
meaout <- cbind(typemarkers,max.df)

meaout$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
meacpm <- ggplot(meaout, aes(y=factor(typemarkers, levels=typemarkers), x=CPM, color=log10(max)))  +theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(1,4.5)) +  theme(axis.line  = element_blank())  +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

plotf <- plot & theme_nothing() & theme(axis.line.x=element_line(color="black", size=1)) & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) & scale_fill_manual(values=c(rep("royalblue3",19), rep("red",15)))


pdf(file="X:/Seurat/MeA_ViolinTest_alternate.pdf", width=60, height=154.8)
(plotf|meacpm) + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()

#pdf(file="X:/Seurat/MeA_ViolinTest_alternatenames.pdf", width=60, height=154.8)
#(plotf|meacpm) + plot_layout(widths=c("1","0.05")) & cowplot::theme_cowplot()
#dev.off()
pdf(file="X:/Seurat/MeA_ViolinTestwide.pdf", width=95, height=37)
(plotf|meacpm) + plot_layout(widths=c("1","0.05")) 
dev.off()

pdf(file="X:/Seurat/MeA_ViolinTest_alternatenames.pdf", width=57, height=154.8)
plotf
dev.off()
pdf(file="X:/Seurat/MeA_AltCPM.pdf", width=3, height=154.8)
meacpm
dev.off()


genes.ph <- CPM[,phenotypes]
genes.ph <- as.matrix(genes.ph)
max <- colMaxs(genes.ph)
max.df <- data.frame(max)
MeAph <- cbind(phenotypes,max.df)
MeAph$CPM <- "CPM"

MeAphCPM <- ggplot(MeAph, aes(y=factor(phenotypes, levels=rev(phenotypes)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(2,3.5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())


plotph <- VlnPlot(MeA, features=phenotypes, pt.size=0, ncol=1)
plotphf <- plotph & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & scale_fill_manual(values=c(rep("royalblue3",19), rep("red",15)))

pdf(file="X:/Seurat/MeA_Phenotypes.pdf", width=50, height=20)
(plotphf|MeAphCPM) + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()




POA <- readRDS("X:/Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(POA)
POA <- RenameIdents(POA, new.cluster.ids)
DefaultAssay(POA) <- "RNA"
POA <- NormalizeData(POA)
POAtypemarkers <- read.table("X:/DotPlotMarkerLists/POA_ReorderedFinalalternate_CPM.txt", header=FALSE)
POAtypemarkers <- unlist(POAtypemarkers)
POAtypemarkers <- unlist(POAtypemarkers)
POAtypemarkers.m <- as.matrix(POAtypemarkers)
POAlevels <- c(1,4,6,10,12,16,18,20,23,25,27,33,39,40,2,3,5,7,8,9,11,13,14,15,17,19,21,22,24,26,28,29,30,31,32,34,35,36,37,38)
#POAlevels <- rev(POAlevels)
Idents(POA) <- factor(Idents(POA), levels=POAlevels)

plot <- VlnPlot(POA, features=rev(POAtypemarkers), pt.size=0, ncol=1)


counts <- POA@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- as.matrix(CPM)
CPM <- t(CPM)

genes.m <- CPM[,POAtypemarkers]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
POAout <- cbind(POAtypemarkers,max.df)

POAout$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
poacpm <- ggplot(POAout, aes(y=factor(POAtypemarkers, levels=POAtypemarkers), x=CPM, color=log10(max)))  +theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(1,4.5)) +  theme(axis.line  = element_blank())  +theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) + theme(axis.title=element_blank()) + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())

plotf <- plot& theme_nothing() & theme(axis.line.x=element_line(color="black", size=1)) & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) &  scale_fill_manual(values=c(rep("royalblue3",14), rep("red",26)))

pdf(file="X:/Seurat/POA_ViolinTest_alternate.pdf", width=60, height=180)
(plotf|poacpm)  + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) 
dev.off()
#pdf(file="X:/Seurat/POA_ViolinTest_alternate.pdf", width=57, height=180)
#plotf & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) 
#dev.off()
pdf(file="X:/Seurat/POA_CPMalt.pdf", width=3, height=180)
poacpm  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA)) 
dev.off()

pdf(file="X:/Seurat/POA_ViolinTest_alternatenames.pdf", width=60, height=180)
(plotf|poacpm)  + plot_layout(widths=c("1","0.05")) & cowplot::theme_cowplot()
dev.off()
pdf(file="X:/Seurat/POA_ViolinTestwide.pdf", width=95, height=43)
(plotf|poacpm)  + plot_layout(widths=c("1","0.05")) 
dev.off()


genes.ph <- CPM[,phenotypes]
genes.ph <- as.matrix(genes.ph)
max <- colMaxs(genes.ph)
max.df <- data.frame(max)
POAph <- cbind(phenotypes,max.df)
POAph$CPM <- "CPM"

POAphCPM <- ggplot(POAph, aes(y=factor(phenotypes, levels=rev(phenotypes)), x=CPM, color=log10(max))) +cowplot::theme_nothing() + geom_point(size=30) + scale_color_viridis(limits=c(2,3.5)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_blank()) + theme(axis.ticks = element_blank())


plotph <- VlnPlot(POA, features=phenotypes, pt.size=0, ncol=1)
plotphf <- plotph & theme(axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks=element_blank(), axis.title=element_blank(), plot.title=element_blank()) & scale_fill_manual(values=c(rep("royalblue3",14), rep("red",26)))

pdf(file="X:/Seurat/POA_Phenotypes.pdf", width=50, height=20)
(plotphf|POAphCPM) + plot_layout(widths=c("1","0.05"))  & theme(panel.background=element_rect(fill="transparent", color=NA), plot.background=element_rect(fill="transparent", color=NA))
dev.off()

