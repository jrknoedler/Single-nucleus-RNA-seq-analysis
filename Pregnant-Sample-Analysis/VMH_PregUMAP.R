#!/usr/bin/env Rscript

output <- "Seurat/VMH_PregIntegrated/VMH_PregFiltered_GrantPlots"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

VMH <- readRDS("Seurat/VMH_PregIntegrated/VMH_PregFiltered_NaiveMergeTest1.rds")
VMH <- subset(VMH, subset=sex=="Male", invert=TRUE)
pdf(file=paste0(output,"UMAP1.pdf"))
DimPlot(VMH, label=TRUE, label.size=9) & theme(legend.position="none")
dev.off()