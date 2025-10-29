#!/usr/bin/env Rscript
output <- "Seurat/VMH_BusMergeTest/VMH_filtered1plots_"

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

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_filtered1.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_Cckar_featurePlot.pdf"))
FeaturePlot(mySeurat, features=c("Cckar"))
dev.off()
pdf(paste0(output, "_Cckar_featurePlot_hormonesplit.pdf"))
FeaturePlot(mySeurat, features=c("Cckar"), split.by="hormone")
dev.off()
pdf(paste0(output, "_Cckar_featurePlot_sexsplit.pdf"))
FeaturePlot(mySeurat, features=c("Cckar"), split.by="sex")
dev.off()
pdf(paste0(output, "_Cyp19a1_featurePlot.pdf"))
FeaturePlot(mySeurat, features=c("Cyp19a1"))
dev.off()
pdf(paste0(output, "_Cyp19a1_featurePlot_sexsplit.pdf"))
FeaturePlot(mySeurat, features=c("Cyp19a1"), split.by="sex")
dev.off()

