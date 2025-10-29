#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MaleVMH_Preprocessing/BNSTPrimed_Hbbtest_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/MaleVMH_Preprocessing/MaleVMHPrimed_Preprocessing.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "Hbbs_1.pdf"))
VlnPlot(mySeurat, features=c("Hbb-bs","Hba-a2","Hbb-ba"), idents=1:7,  pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "Hbbs_2.pdf"))
VlnPlot(mySeurat, features=c("Hbb-bs","Hba-a2","Hbb-ba"), idents=8:14, ,pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "Hbbs_3.pdf"))
VlnPlot(mySeurat, features=c("Hbb-bs","Hba-a2","Hbb-ba"), idents=15:21,pt.size=0, ncol=1)
dev.off()
pdf(paste0(output, "Hbbs_4.pdf"))
VlnPlot(mySeurat, features=c("Hbb-bs","Hbb-bt","Hbb-ba"), idents=22:28,pt.size=0, ncol=1)
dev.off()

