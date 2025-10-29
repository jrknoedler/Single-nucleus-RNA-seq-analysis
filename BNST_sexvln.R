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

output <- "Seurat/BNST_IndependentAnalysis/BNST_Figs5heatmapsfinal_"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
Sexmarkers <- read.table("DotPlotMarkerLists/BNST_Sexmarkers.txt")
Sexmarkers <- unlist(Sexmarkers)
Sexclusters <- read.table("DotPlotMarkerLists/BNST_Markedclusters.txt")
Sexclusters <- unlist(Sexclusters)
pdf(file=paste0(output,"BNST_sexmarkerVlns.pdf"), width=20,height=100)
VlnPlot(mySeurat, features=Sexmarkers, idents=Sexclusters, ncol=1, pt.size=0)
dev.off()