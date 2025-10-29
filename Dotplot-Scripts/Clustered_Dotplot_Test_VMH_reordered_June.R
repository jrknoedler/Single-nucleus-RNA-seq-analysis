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

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
a <- AverageExpression(mySeurat, assays=c("RNA"), features=c("Esr1", "Ar","Pgr","Esr2"))
write.csv(a, file=paste0(output,"Esr1perclust.csv"))
typemarkers <- read.table("DotPlotMarkerLists/VMH_reorderedfinal.txt", header=FALSE)
typemarkers2 <- read.table("DotPlotMarkerLists/VMH_reorderedfinal_CPM.txt", header=FALSE)
typemarkers2 <- unlist(typemarkers2)
typemarkers <- unlist(typemarkers)
typemarkers.m <- as.matrix(typemarkers)
mySeurat <- ScaleData(mySeurat)


mylevels <- c(2,3,4,6,7,11,12,15,16,17,19,20,21,24,25,26,27,1,5,8,9,10,13,14,18,22,23)
mylevels <- rev(mylevels)
pdf(file=paste0(output,"markerVln.pdf"), width=30, height=20)
VlnPlot(mySeurat, features=c("Esr1","Pgr","Slc17a6","Cckar"), pt.size=0, ncol=1)
dev.off()

pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers.m, dot.min=0.25)
dev.off()

Pct <- DotPlot(mySeurat, features=typemarkers.m)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
pdf(file=paste0(output,"_rawdotplotscaled.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"))+ scale_size(limits=c(25,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()
genedots <- data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) +  cowplot::theme_nothing() +
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"))+ scale_size(limits=c(25,100)) +
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_blank(), axis.text.y=element_text()) +
  ylab('') +
  theme(axis.ticks = element_blank())

counts <- mySeurat@assays[["RNA"]]@counts
counts.m <- as.matrix(counts)
CPM <- RelativeCounts(counts, scale.factor=1e6, verbose=TRUE)
CPM <- t(CPM)

genes.m <- CPM[,typemarkers2]
dim(genes.m)
head(genes.m)
genes.f <- as.matrix(genes.m)
max <- colMaxs(genes.f)
max.df <- data.frame(max)
out <- cbind(typemarkers2,max.df)
head(out)
out$CPM <- "CPM"
#typemarkers <- typemarkers[,data$features.plot[["order"]]]
write.csv(out, file=paste0(output,"CPMcompiled.csv"))
p <- ggplot(out, aes(x=factor(typemarkers2, levels=typemarkers2), y=CPM, color=max)) +cowplot::theme_nothing() + geom_point(size=5) + scale_color_viridis() +  theme(axis.line  = element_blank())   + theme(axis.title=element_blank()) + labs(x="") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(axis.ticks = element_blank())
pdf(paste0(output,"Cpmdotdemo.pdf"), width=16)
p
dev.off()
pdf(paste0(output,"combinedottest.pdf"), width=16)
(genedots/p) + plot_layout(heights=c("1","0.05"))
dev.off()
