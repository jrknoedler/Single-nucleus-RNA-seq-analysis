#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_sharedmarkers"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

filteredSeurat <- readRDS("Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_doubletsremoved_filtered2.rds")
DefaultAssay(filteredSeurat) <- "RNA"
NormalizeData(filteredSeurat)
for (i in 0:42){
try({
filteredSeurat.markers <- FindConservedMarkers(filteredSeurat, ident.1= i, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, grouping.var = "sex")
write.csv(filteredSeurat.markers, file=paste0(output,i,"_allposmarkers.csv"))
})
}

