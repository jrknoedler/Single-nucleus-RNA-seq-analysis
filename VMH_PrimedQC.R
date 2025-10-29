#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMHprimed_exampleQC"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentRound1__Primed.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)

pdf(file=paste0(output,"Esr1ftr.pdf"))
FeaturePlot(mySeurat, features=c("Esr1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output,"Sypftr.pdf"))
FeaturePlot(mySeurat, features=c("Syp"), order=TRUE, cols=c("light blue","red"))
dev.off()