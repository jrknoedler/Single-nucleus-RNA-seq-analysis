#!/usr/bin/env Rscript

output <- "Seurat/VMH_Anderson/FilteredQC_"

library(Seurat)
#library(cowplot)
#library(reticulate)
#library(ggplot2)
library(dplyr)
#library(MAST)
library(future)
library(Matrix)

mySeurat <- readRDS("Seurat/VMH_Anderson/Anderson_NonneuronalFilt5.rds")
pdf(file=paste0(output,"UMI.pdf"), width=100, height=40)
VlnPlot(mySeurat, features=c("nCount_RNA"), pt.size=0, split.by="Condition")
dev.off()
pdf(file=paste0(output,"Genes.pdf"), width=100, height=40)
VlnPlot(mySeurat, features=c("nFeature_RNA"), pt.size=0, split.by="Condition")
dev.off()