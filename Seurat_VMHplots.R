#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_MalevPrimedMerged_Exploratory/VMH_Miscgenes_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/VMH_MalevPrimedMerged_Exploratory/VMH_MalevPrimed_Exploratory_filtered2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_Vln_Markers1.pdf"))
VlnPlot(mySeurat, features =c("Esr1", "Satb2","Dlk1"),ncol=1, split.by="sex", pt.size=0)
dev.off()
pdf(paste0(output, "_Vln_Markers2.pdf"))
VlnPlot(mySeurat, features =c("Pgr", "Setdb2","Cckar"),ncol=1, split.by="sex", pt.size=0)
dev.off()
pdf(paste0(output, "_Vln_Markers3.pdf"))
VlnPlot(mySeurat, features =c("Tac1", "Avp","Cck"),ncol=1, split.by="sex", pt.size=0)
dev.off()
pdf(paste0(output, "_Vln_Markers4.pdf"))
VlnPlot(mySeurat, features =c("Rprm", "Pdyn","Cyp19a1"),ncol=1, split.by="sex", pt.size=0)
dev.off()
