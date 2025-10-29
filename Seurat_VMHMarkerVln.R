#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_BusMergeTest/MarkerVln/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_ExcludeDE.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
pdf(file=paste0(output,"_Slc17a6.pdf"), width=20, height=10)
vlnvar <- VlnPlot(mySeurat, features="Slc17a6", pt.size=0, combine=TRUE)
plot(vlnvar)
dev.off()
genelist <- read.table("Genelists/VMHMarkers.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
for (g in genelist)
try({
pdf(file=paste0(output,g,".pdf"), width=20, height=10)
vlnvar <- VlnPlot(mySeurat, features=c(g),, pt.size=0, combine=TRUE)
plot(vlnvar)
graphics.off()
})



