#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNSTUnprimed_exampleQC"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentRound1__Unprimed.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)

pdf(file=paste0(output,"Esr1ftr.pdf"))
FeaturePlot(mySeurat, features=c("Esr1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output,"Sypftr.pdf"))
FeaturePlot(mySeurat, features=c("Syp"), order=TRUE, cols=c("light blue","red"))
dev.off()