#!/usr/bin/env Rscript

output <- "Seurat/UnprimedMeA_Exploratory/UnprimedMeA_1000cutoff_filtered3"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/UnprimedMeA_Exploratory/UnprimedMeA_1000cutoff_filtered3.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output,"_Esr1Vln.pdf"))
VlnPlot(mySeurat, features="Esr1", pt.size=0)
dev.off()
pdf(paste0(output,"_Slc17a6Vln.pdf"))
VlnPlot(mySeurat, features="Slc17a6", pt.size=0)
dev.off()
pdf(paste0(output,"_Gad2Vln.pdf"))
VlnPlot(mySeurat, features="Gad2", pt.size=0)
dev.off()
pdf(paste0(output,"_Slc32a1Vln.pdf"))
VlnPlot(mySeurat, features="Slc32a1", pt.size=0)
dev.off()
pdf(paste0(output,"_UMIVln.pdf"))
VlnPlot(mySeurat, features="nCount_RNA", pt.size=0)
dev.off()



