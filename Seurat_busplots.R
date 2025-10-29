#!/usr/bin/env Rscript
output <- "Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_excludeDE"

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

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_ExcludeDE.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_Cckar_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Cckar"))
dev.off()

