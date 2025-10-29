#!/usr/bin/env Rscript
library(Seurat)
library(rgeos)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)
library(SeuratDisk)


mySeurat <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")
mySeurat.loom <- as.loom(mySeurat, filename="loom/POA_PregMerged_Updated.loom")
mySeurat.loom$close_all()