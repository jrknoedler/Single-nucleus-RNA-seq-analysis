#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_Plots_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")
pdf(file=paste0(output,"CartptVln.pdf"), width=14)
VlnPlot(mySeurat, features=c("Cartpt", "Tac1"), ncol=1, pt.size=0)
dev.off()
pdf(file=paste0(output,"CartptFeaturePlot.pdf"), width=14)
FeaturePlot(mySeurat, features=c("Cartpt"), cols=c("light blue", "red"))
dev.off()
