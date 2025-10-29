#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_MalevPrimedMerged_Exploratory/VMH_CCAtest_ConservedMarkers_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/VMH_MalevPrimedMerged_Exploratory/VMH_MalevPrimed_Exploratory_filtered1.rds")
mySeurat
DefaultAssay(mySeurat) <- "RNA"
NormalizeData(mySeurat)
for (i in 0:34){
try({
mySeurat.markers <- FindConservedMarkers(mySeurat, ident.1= i, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, grouping.var = "sex")
write.csv(mySeurat.markers, file=paste0(output,i,"_allposmarkers.csv"))
})
}
saveRDS(mySeurat, file=(paste0(output, ".rds")))
