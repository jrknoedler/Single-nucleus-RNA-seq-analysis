#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/MeA_Celltypes_reorderedJune"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
a <- AverageExpression(mySeurat, assays=c("RNA"), features=c("Esr1", "Ar","Pgr","Esr2"))
write.csv(a, file=paste0(output,"Esr1perclust.csv"))
typemarkers <- read.table("DotPlotMarkerLists/MeA_ReorderedFinal.txt", header=FALSE)
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
mySeurat <- ScaleData(mySeurat)


mylevels <- c(1,4,9,13,14,16,18,21,22,23,24,29,30,31,32,33,2,3,5,6,8,11,12,17,19,20,25,26,27,28,34,7,10,15)
mylevels <- rev(mylevels)
pdf(file=paste0(output,"markerVln.pdf"), width=30, height=20)
VlnPlot(mySeurat, features=c("Esr1","Pgr","Slc17a6","Cckar"), pt.size=0, ncol=1)
dev.off()

pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.25)
dev.off()

Pct <- DotPlot(mySeurat, features=typemarkers)
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

