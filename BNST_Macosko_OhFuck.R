#!/usr/bin/env Rscript

output <- "Seurat/BNSTMacosko/Macosko_SanityCheck"
library(Seurat)
#library(cowplot)
#library(reticulate)
#library(ggplot2)
library(dplyr)
#library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)

macosko <- readRDS("Seurat/BNSTMacosko/Macosko_BNST_Exploratory_res1_annotated.rds")
DefaultAssay(macosko) <- "RNA"
macosko <- NormalizeData(macosko)
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=45)
VlnPlot(macosko, features=c("Gad1","Slc17a6","Esr1","Pgr","Ar","Prlr","Cyp19a1","Tac1"), ncol=1, pt.size=0)
dev.off()