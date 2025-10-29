#!/usr/bin/env Rscript


output <- "Seurat/BNST_IndependentAnalysis/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered.rds")
mat <- mySeurat[["SCT"]]@scale.data

mat <- t(mat)

write.csv(mat, file=paste0(output,"vst_countstrans.csv"))
