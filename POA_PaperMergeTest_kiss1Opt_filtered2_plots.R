#!/usr/bin/env Rscript

output <- "Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filteredplots"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered2.rds")

pdf(file=paste0(output,"Kiss1_Ftplot.pdf"))
FeaturePlot(mySeurat, features=c("Kiss1"))
dev.off()

pdf(file=paste0(output,"Kiss1_Ftplot_hormonesplit.pdf"), width=12)
FeaturePlot(mySeurat, features=c("Kiss1"), split.by="Hormone")
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))