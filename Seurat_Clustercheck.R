#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_Doubletsremoved/MalevUnprimed_Fineclusters_sharedmarkers_filtered_doubletsremoved"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
Cluster.diff1 <- FindMarkers(mySeurat, ident.1 = "25", ident.2="38", min.pct=0.1, logfc.threshold=0.25, verbose=TRUE)
write.csv(Cluster.diff1, file=paste0(output, "_cluster25v38.csv"))
Cluster.diff2 <- FindMarkers(mySeurat, ident.1 = "6", ident.2="26", min.pct=0.1, logfc.threshold=0.25, verbose=TRUE)
write.csv(Cluster.diff2, file=paste0(output, "_cluster6v26.csv"))
Cluster.diff3 <- FindMarkers(mySeurat, ident.1 = "2", ident.2="37", min.pct=0.1, logfc.threshold=0.25, verbose=TRUE)
write.csv(Cluster.diff3, file=paste0(output, "_cluster2v37.csv"))
