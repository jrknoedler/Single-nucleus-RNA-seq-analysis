#!/usr/bin/env Rscript

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
output <- "Seurat/POA_MalevUnprimedMerged_Exploratory/Vln"

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_Chungha1.pdf"))
VlnPlot(mySeurat, features=c("Drd1","Drd2","Htr2a"), pt.size=0, idents=c(0:9), split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "_Chungha2.pdf"))
VlnPlot(mySeurat, features=c("Drd1","Drd2","Htr2a"), pt.size=0, idents=c(10:19), split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "_Chungha3.pdf"))
VlnPlot(mySeurat, features=c("Drd1","Drd2","Htr2a"), pt.size=0, idents=c(20:29), split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "_Chungha4.pdf"))
VlnPlot(mySeurat, features=c("Drd1","Drd2","Htr2a"), pt.size=0, idents=c(30:38), split.by="sex", ncol=1)
dev.off()
