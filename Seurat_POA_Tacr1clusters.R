#!/usr/bin/env Rscript
output <- "Seurat/POA_BusMergeTest/POA_Busplotstacr1clusters_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
library(matrixStats)
library(qlcMatrix)
library(igraph)
library(RANN)

mySeurat <- readRDS("Seurat/POA_BusMergeTest/POA_BusMergeTest_tacr1filtered.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_Tacr1_featurePlot.pdf"))
VlnPlot(mySeurat, features=c("Tacr1"), ncol=1, pt.size=0)
dev.off()


