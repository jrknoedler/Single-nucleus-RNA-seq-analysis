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

output <- "Seurat/POA_IndependentAnalysis/POA_Figs5heatmapsfinal_"

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
Sexmarkers <- read.table("DotPlotMarkerLists/POA_Sexmarkers.txt")
Sexmarkers <- unlist(Sexmarkers)
Sexclusters <- read.table("DotPlotMarkerLists/POA_Markedclusters.txt")
Sexclusters <- unlist(Sexclusters)
pdf(file=paste0(output,"POA_sexmarkerVlns.pdf"), width=20,height=100)
VlnPlot(mySeurat, features=Sexmarkers, idents=Sexclusters, ncol=1, pt.size=0)
dev.off()