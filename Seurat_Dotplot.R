
#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_Doubletsremoved/DiffgenesClustExp/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
idents <- c(0:9)
subset < - subset(mySeurat, idents=idents)
markers.to.plot <- c("Esr1","Pappa","Kirrel","Pdzrn4","Ar")
pdf(paste0(output, "dotplot.pdf"))
DotPlot(subset, features=markers.to.plot, split.by="sex", cols = c("blue","red")) + RotatedAxis()
dev.off()
