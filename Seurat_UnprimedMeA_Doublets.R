#!/usr/bin/env Rscript

output <- "Seurat/UnprimedMeA_Exploratory/UnprimedMeA_emptyDrops_doublets_30pcs_lowfdr"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DoubletFinder)

mySeurat <- readRDS("Seurat/UnprimedMeA_Exploratory/UnprimedMeA_emptyDrops_doublets_30pcs_lowfdr.rds")
pdf(file=paste0(output,"doublettest.pdf"))
DimPlot(mySeurat, group.by="DF.classifications_0.25_0.14_1337")
dev.off()