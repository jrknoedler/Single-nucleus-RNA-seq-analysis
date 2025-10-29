#!/usr/bin/env Rscript

output <- "Seurat/VMH_BusMergeTest/VMH_featureplots"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_filtered1.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output,"_Xist.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Xist", split.by="status")
dev.off()
pdf(paste0(output,"_Uty.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Uty", split.by="status")
dev.off()

