#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/FullPOA_Integrated/FullPOA_finalfiltered_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(clustree)

mySeurat <- readRDS("Seurat/FullPOA_Integrated/FullPOA_furtherfilteredunnormalized.rds")
mySeurat
mySeurat <- BuildClusterTree(mySeurat, reorder=TRUE, verbose=TRUE)
pdf(paste0(output, "clustree.pdf"))
PlotClusterTree(object=mySeurat)
dev.off()

