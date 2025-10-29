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

output <- "Seurat/MeA_IndependentAnalysis/MeA_Figs5heatmapsfinal_"

mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
Sexmarkers <- read.table("DotPlotMarkerLists/MeA_Sexmarkers.txt")
Sexmarkers <- unlist(Sexmarkers)
Sexclusters <- read.table("DotPlotMarkerLists/MeA_Markedclusters.txt")
Sexclusters <- unlist(Sexclusters)
pdf(file=paste0(output,"MeA_sexmarkerVlns.pdf"), width=20,height=100)
VlnPlot(mySeurat, features=Sexmarkers, idents=Sexclusters, ncol=1, pt.size=0)
dev.off()