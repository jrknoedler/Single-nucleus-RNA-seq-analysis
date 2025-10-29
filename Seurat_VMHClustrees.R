#!/usr/bin/env Rscript

output <- "Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_clustree_genes_"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(clustree)

mySeurat <- readRDS("Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_filtered1_clustree.rds")
pdf(file=paste0(output,"_tree_Esr1.pdf"))
clustree(mySeurat, prefix="SCT_snn_res.", node_colour="Esr1", node_colour_aggr="max")
dev.off()
pdf(file=paste0(output,"_tree_Slc17a6.pdf"))
clustree(mySeurat, prefix="SCT_snn_res.", node_colour="Slc16a6", node_colour_aggr="max")
dev.off()
pdf(file=paste0(output,"_tree_Gad1.pdf"))
clustree(mySeurat, prefix="SCT_snn_res.", node_colour="Gad1", node_colour_aggr="max")
dev.off()
