#!/usr/bin/env Rscript

output <- "Seurat/Pregnancy_Full_Dataset_Analysis/Merged_Regulon_Analysis"

library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
#library(SCENIC)
#library(Matrix)
#library(AUCell)
#library(loomR)
library(future)

options(future.globals.maxSize=8000*1024^2)
mySeurat <- readRDS("Seurat/Pregnancy_Full_Dataset_Analysis/Merged_Pregnancy_with_Positive_Regulons_1stpass_regsonly.rds")
head(mySeurat[[]])

DefaultAssay(mySeurat) <- "Regulons"
mySeurat <- SetIdent(mySeurat, value=mySeurat@meta.data$Region)
head(Idents(mySeurat))
region_reg <- FindAllMarkers(mySeurat, test.use="roc", verbose=TRUE, logfc.threshold=0, , return.thresh=0)
write.csv(region_reg, file=paste0(output,"Regulon_ROC_Regionvar.csv"))
