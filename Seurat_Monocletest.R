#!/usr/bin/env Rscript
output <- "Seurat/BNST_BusMergeTest/BNST_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
library(matrixStats)
library(qlcMatrix)
library(igraph)
library(RANN)

mySeurat <- readRDS("Seurat/BNST_BusMergeTest/BNST_BusMergeTest_nocastrate_filtered1.rds")
monocle <- importCDS(mySeurat)
monocle

