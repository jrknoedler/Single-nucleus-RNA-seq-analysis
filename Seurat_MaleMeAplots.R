#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MaleMeA_Preprocessing/Aromataseplots"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/MaleMeA_Preprocessing/MaleMeA_Preprocessing.rds")
mySeurat <- subset(mySeurat, idents=c("12","13","16","21","24","27"))
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_Vln_Markers1.pdf"))
VlnPlot(mySeurat, features =c("Esr2", "Cyp19a1","Ar"),ncol=1, pt.size=0)
dev.off()
pdf(paste0(output, "_Vln_Markers2.pdf"))
VlnPlot(mySeurat, features =c("Prlr", "Tacr1","Cckar"),ncol=1, pt.size=0)
dev.off()
pdf(paste0(output, "_Vln_Markers3.pdf"))
VlnPlot(mySeurat, features =c("Tac1", "Prl","Cck"),ncol=1, pt.size=0)
dev.off()
pdf(paste0(output, "_Vln_Markers4.pdf"))
VlnPlot(mySeurat, features =c("St18", "Moxd1","Pappa"),ncol=1, pt.size=0)
dev.off()
