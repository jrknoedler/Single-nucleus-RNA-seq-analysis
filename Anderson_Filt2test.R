#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_Anderson/Anderson_NonneuronalFilt2_markers"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)


mySeurat <- readRDS("Seurat/VMH_Anderson/Anderson_NonneuronalFilt2.rds")
mySeurat <- NormalizeData(mySeurat)
pdf(file=paste0(output,".pdf"), width=50, height=12)
VlnPlot(mySeurat, features=c("Slc17a6","Slc32a1","Esr1","Cckar"), ncol=1, pt.size=0)
dev.off()

