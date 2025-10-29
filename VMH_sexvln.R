#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"


library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)

output <- "Seurat/VMH_IndependentAnalysis/VMH_Figs5heatmapsfinal_"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
Sexmarkers <- read.table("DotPlotMarkerLists/VMH_Sexmarkers.txt")
Sexmarkers <- unlist(Sexmarkers)
Sexclusters <- read.table("DotPlotMarkerLists/VMH_Markedclusters.txt")
Sexclusters <- unlist(Sexclusters)
pdf(file=paste0(output,"VMH_sexmarkerVlns.pdf"), width=20,height=100)
VlnPlot(mySeurat, features=Sexmarkers, idents=Sexclusters, ncol=1, pt.size=0)
dev.off()
pdf(file=paste0(output,"VMH_sexmarkerVlnsall.pdf"), width=30,height=100)
VlnPlot(mySeurat, features=Sexmarkers, ncol=1, pt.size=0)
dev.off()