#!/usr/bin/env Rscript

output <- "Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_excludeDE_clustree"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(clustree)

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_ExcludeDE.rds")
mySeurat <- FindClusters(mySeurat, resolution=c(0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0), plot.SNN=TRUE, save.SNN=TRUE)
head(mySeurat[[]])
pdf(file=paste0(output,"_tree1.pdf"))
clustree(mySeurat, prefix="SCT_snn_res.")
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))
