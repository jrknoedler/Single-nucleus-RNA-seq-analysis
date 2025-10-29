#!/usr/bin/env Rscript

output <- "Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_excludeDE_MAST"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(clustree)
library(MAST)

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_ExcludeDE.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
head(mySeurat[[]])
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))

