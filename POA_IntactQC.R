#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POAMale_exampleQC"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Intact_Kiss1opt_500_Intact.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)

pdf(file=paste0(output,"Esr1ftr.pdf"))
FeaturePlot(mySeurat, features=c("Esr1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output,"Sypftr.pdf"))
FeaturePlot(mySeurat, features=c("Syp"), order=TRUE, cols=c("light blue","red"))
dev.off()