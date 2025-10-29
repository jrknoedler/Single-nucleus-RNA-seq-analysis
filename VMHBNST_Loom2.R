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


#mySeurat <- readRDS("Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")
#mySeurat.loom <- as.loom(mySeurat, filename="loom/VMH_PregMerged_Updated.loom")
#mySeurat.loom$close_all()
mySeurat2 <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered1.rds")
mySeurat2.loom <- as.loom(mySeurat2, filename="loom/BNST_PregMerged_Updated.loom")
mySeurat2.loom$close_all()