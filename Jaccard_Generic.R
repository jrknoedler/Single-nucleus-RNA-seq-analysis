#!/usr/bin/env Rscript
rm(list=ls())

library(scclusteval)
library(Seurat)
library(tidyverse)
library(dplyr)

library(future)
library(Matrix)

options(future.globals.maxSize=8000*1024^2)


### VMH ====

output <- "/scratch/users/knoedler/Seurat/MeA_Barcode_Transfer/Jaccard_Plots"

##Load dataset 1 (Published data)
Published <- readRDS("/scratch/users/knoedler/Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
Published <- RenameCells(object = Published, add.cell.id = "Ref")
Published@meta.data$NewIntegration <- "Published"
Published@meta.data$PublishedCluster <- Idents(Published)

##Load dataset2 (Pregnancy Integration)
Integration <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/VMH/VMH_Naive_Integration_PC30_Res12_Fixed.rds")
Integration@meta.data$IntegrationCluster <- Idents(Integration)

##Jaccard Index

#Jaccard Similarity Heatmap
pdf(file=paste0(output,"_PairWiseJaccardSetsHeatmap.pdf"), width=27, height=25)
PairWiseJaccardSetsHeatmap(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                           set_names(Integration@meta.data$IntegrationCluster, nm=colnames(Integration)),
                           show_row_dend = F, show_column_dend = F, cluster_row = F, cluster_column =F,
                           col_low = "#FFFFFF", col_high = "#004D40")
dev.off()

#Jaccard Index Calculation
PairWiseJaccardSets <- PairWiseJaccardSets(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                                           set_names(Integration@meta.data$IntegrationCluster, nm=colnames(Integration)))
write.csv(PairWiseJaccardSets, file = paste0(output, "_PairWiseJaccardSets.csv"))
