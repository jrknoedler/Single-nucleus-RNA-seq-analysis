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

VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(VMH)
VMH <- RenameIdents(VMH, new.cluster.ids)
DefaultAssay(VMH) <- "RNA"
VMH <- NormalizeData(VMH)
a <- AverageExpression(VMH, assays=c("RNA"), features=c("Esr1", "Ar","Pgr","Esr2"))
write.csv(a, file=paste0(output,"Esr1perclust.csv"))
vmhtypemarkers <- read.table("DotPlotMarkerLists/VMH_reorderedfinal.txt", header=FALSE)
vmhtypemarkers2 <- read.table("DotPlotMarkerLists/VMH_reorderedfinal_CPM.txt", header=FALSE)
vmhtypemarkers2 <- unlist(vmhtypemarkers2)
vmhtypemarkers <- unlist(vmhtypemarkers)
vmhtypemarkers.m <- as.matrix(vmhtypemarkers)
MeA <- ScaleData(VMH)


VMHlevels <- c(2,3,4,6,7,11,12,15,16,17,19,20,21,24,25,26,27,1,5,8,9,10,13,14,18,22,23)
VMHlevels <- rev(VMHlevels)



Pct <- DotPlot(VMH, features=vmhtypemarkers.m)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

VMHdata <- Pct$data


VMHdots <- VMHdata %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=VMHlevels), color=avg.exp.scaled, size=pct.exp)) +  cowplot::theme_cowplot() +
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5,2.5))+ scale_size(limits=c(25,100),breaks=seq(25,100,25)) +
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

counts <- VMH@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,vmhtypemarkers2]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
vmhout <- cbind(vmhtypemarkers2,max.df)
write.csv(vmhout, file="Seurat/VMH_Markercpm.csv")
vmhout$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
vmhcpm <- ggplot(vmhout, aes(x=factor(vmhtypemarkers2, levels=vmhtypemarkers2), y=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=5) + scale_color_viridis(limits=c(2.2,4.2)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.ticks = element_blank())


VMHplot <- (VMHdots/vmhcpm) + plot_layout(heights=c("1","0.05"))



MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
names(new.cluster.ids) <- levels(MeA)
MeA <- RenameIdents(MeA, new.cluster.ids)
DefaultAssay(MeA) <- "RNA"
MeA <- NormalizeData(MeA)


typemarkers <- read.table("DotPlotMarkerLists/MeA_ReorderedFinal.txt", header=FALSE)
typemarkers2 <- read.table("DotPlotMarkerLists/MeA_ReorderedFinal_CPM.txt", header=FALSE)
typemarkers2 <- unlist(typemarkers2)
typemarkers <- unlist(typemarkers)
typemarkers.m <- as.matrix(typemarkers)

MeA <- ScaleData(MeA)


MeAlevels <- c(1,4,9,13,14,16,18,21,22,23,24,29,30,31,32,33,2,3,5,6,8,11,12,17,19,20,25,26,27,28,34,7,10,15)
MeAlevels <- rev(MeAlevels)




Pct <- DotPlot(MeA, features=typemarkers.m)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

MeAdata <- Pct$data


MeAdots <- MeAdata %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=MeAlevels), color=avg.exp.scaled, size=pct.exp)) +  cowplot::theme_cowplot() +
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5,2.5))+ scale_size(limits=c(25,100)) +
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

counts <- MeA@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,typemarkers2]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
meaout <- cbind(typemarkers2,max.df)
write.csv(meaout, file="Seurat/MeA_Markercpm.csv")
meaout$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
meacpm <- ggplot(meaout, aes(x=factor(typemarkers2, levels=typemarkers2), y=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=5) + scale_color_viridis(limits=c(2.2,4.2)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.ticks = element_blank())


MeAplot <- (MeAdots/meacpm) + plot_layout(heights=c("1","0.05"))
pdf("Seurat/Fig4test.pdf", width=40)
(VMHplot|MeAplot) + plot_layout(widths=c("0.7","1"))
dev.off()


POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(POA)
POA <- RenameIdents(POA, new.cluster.ids)
DefaultAssay(POA) <- "RNA"
POA <- NormalizeData(POA)


POAtypemarkers <- read.table("DotPlotMarkerLists/POA_ReorderedFinal.txt", header=FALSE)
POAtypemarkers2 <- read.table("DotPlotMarkerLists/POA_ReorderedFinal_CPM.txt", header=FALSE)
POAtypemarkers2 <- unlist(POAtypemarkers2)
POAtypemarkers <- unlist(POAtypemarkers)
POAtypemarkers.m <- as.matrix(POAtypemarkers)

POA <- ScaleData(POA)


POAlevels <- c(1,4,6,10,12,16,18,20,23,25,27,33,39,40,2,3,5,7,8,9,11,13,14,15,17,19,21,22,24,26,28,29,30,31,32,34,35,36,37,38)
POAlevels <- rev(POAlevels)




Pct <- DotPlot(POA, features=POAtypemarkers.m)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data


POAdata <- Pct$data

POAdots <- POAdata %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=POAlevels), color=avg.exp.scaled, size=pct.exp)) +  cowplot::theme_cowplot() +
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5,2.5))+ scale_size(limits=c(25,100)) +
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

counts <- POA@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,POAtypemarkers2]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
POAout <- cbind(POAtypemarkers2,max.df)
write.csv(POAout, file="Seurat/POA_Markercpm.csv")
POAout$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
poacpm <- ggplot(POAout, aes(x=factor(POAtypemarkers2, levels=POAtypemarkers2), y=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=5) + scale_color_viridis(limits=c(2.2,4.2)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.ticks = element_blank())


POAplot <- (POAdots/poacpm) + plot_layout(heights=c("1","0.05"))




BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(BNST)
BNST <- RenameIdents(BNST, new.cluster.ids)
DefaultAssay(BNST) <- "RNA"
BNST <- NormalizeData(BNST)


BNSTtypemarkers <- read.table("DotPlotMarkerLists/BNST_ReorderedFinal.txt", header=FALSE)
BNSTtypemarkers2 <- read.table("DotPlotMarkerLists/BNST_ReorderedFinal_CPM.txt", header=FALSE)
BNSTtypemarkers2 <- unlist(BNSTtypemarkers2)
BNSTtypemarkers <- unlist(BNSTtypemarkers)
BNSTtypemarkers.m <- as.matrix(BNSTtypemarkers)

BNST <- ScaleData(BNST)


BNSTlevels <- c(6,15,20,36,1,2,3,4,5,7,8,9,10,11,12,13,14,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)
BNSTlevels <- rev(BNSTlevels)




Pct <- DotPlot(BNST, features=BNSTtypemarkers.m)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data


BNSTdata <- Pct$data

BNSTdots <- BNSTdata %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=BNSTlevels), color=avg.exp.scaled, size=pct.exp)) +  cowplot::theme_cowplot() +
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5,2.5))+ scale_size(limits=c(25,100)) +
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text(hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())

counts <- BNST@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,BNSTtypemarkers2]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
BNSTout <- cbind(BNSTtypemarkers2,max.df)
write.csv(BNSTout, file="Seurat/BNST_Markercpm.csv")
BNSTout$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
BNSTcpm <- ggplot(BNSTout, aes(x=factor(BNSTtypemarkers2, levels=BNSTtypemarkers2), y=CPM, color=log10(max))) +cowplot::theme_cowplot() + geom_point(size=5) + scale_color_viridis(limits=c(2.2,4.2)) +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.ticks = element_blank())

BNSTplot <- (BNSTdots/BNSTcpm) + plot_layout(heights=c("1","0.05"))

pdf(file="Seurat/Fullfig4_checkscales.pdf", width=82, height=8)
(BNSTplot|MeAplot|POAplot|VMHplot)+plot_layout(widths=c("1","0.94","1.1","0.75"))
dev.off()
