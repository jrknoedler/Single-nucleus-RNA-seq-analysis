#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_Cyp19a1_moremarkers"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_Aromataseneurons.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_set1.pdf"))
VlnPlot(mySeurat, features=c("Nfix","Nfib","Nfia"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "set2.pdf"))
VlnPlot(mySeurat, features=c("Haus4","Fstl4","Kif6"),pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "set3.pdf"))
VlnPlot(mySeurat, features=c("Arhgap36","C79798","Stk32a"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "set4.pdf"))
VlnPlot(mySeurat, features=c("Chd7","Sytl5","Htr2c"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "set5.pdf"))
VlnPlot(mySeurat, features=c("Adamts19","St8sia6","Phactr2"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "set6.pdf"))
VlnPlot(mySeurat, features=c("Trhr","Prox1","Vav3"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "set7.pdf"))
VlnPlot(mySeurat, features=c("Greb1","Syne2","Il1rap"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "set8.pdf"))
VlnPlot(mySeurat, features=c("Sp9","Olfm3","Myo16"), pt.size=0, ncol=1)
dev.off()
