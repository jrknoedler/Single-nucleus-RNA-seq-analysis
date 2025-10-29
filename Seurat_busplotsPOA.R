#!/usr/bin/env Rscript
output <- "Seurat/POA_BusMergeTest/POA_"

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

mySeurat <- readRDS("Seurat/POA_BusMergeTest/POA_BusMergeTest_filtered1.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_Vln_Cyp19_Tacr1.pdf"))
VlnPlot(mySeurat, features="Cyp19a1", pt.size=1, idents="20", split.by="status")
dev.off()
pdf(paste0(output,"Cyp19a1_featureplot.pdf"))
FeaturePlot(mySeurat, features="Cyp19a1")
dev.off()
pdf(paste0(output,"Tacr1_featureplot.pdf"))
FeaturePlot(mySeurat, features="Tacr1")
dev.off()
