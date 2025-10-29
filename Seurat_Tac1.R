#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_Brs3"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_doubletsremoved_filtered2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "Genes_1.pdf"))
VlnPlot(mySeurat, features=c("Brs3","Esr1","Cyp19a1"), idents=0:10,  pt.size=0, ncol=1, split.by="sex")
dev.off()
pdf(paste0(output, "Genes_2.pdf"))
VlnPlot(mySeurat, features=c("Brs3","Esr1","Cyp19a1"), idents=11:20, ,pt.size=0, ncol=1, split.by="sex")
dev.off()
pdf(paste0(output, "Genes_3.pdf"))
VlnPlot(mySeurat, features=c("Brs3","Esr1","Cyp19a1"), idents=21:30,pt.size=0, ncol=1, split.by="sex")
dev.off()
pdf(paste0(output, "Genes_4.pdf"))
VlnPlot(mySeurat, features=c("Brs3","Esr1","Cyp19a1"), idents=31:40,pt.size=0, ncol=1, split.by="sex")
dev.off()

