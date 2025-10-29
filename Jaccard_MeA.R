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
PubIDs <- read.table("Seurat/MeA_Barcode_Transfer/Pub_IDs.txt", header=FALSE)
PubIDs <- unlist(PubIDs)
Published <- readRDS("/scratch/users/knoedler/Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
Published <- RenameCells(object = Published, new.names = PubIDs)
Published@meta.data$NewIntegration <- "Published"
Published@meta.data$PublishedCluster <- Idents(Published)


##Load dataset2 (Pregnancy Integration)
Integration1 <- readRDS("/scratch/users/knoedler/Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_Naive.rds")
Integration1@meta.data$IntegrationCluster <- Idents(Integration1)

##Load dataset3 (Pregnancy Integrationv2)
SCTIDs <- read.table("Seurat/MeA_Barcode_Transfer/SCT_IDs.txt", header=FALSE)
SCTIDs <- unlist(SCTIDs)
Integration2 <- readRDS("/scratch/users/knoedler/Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_SCT.rds")
Integration2 <- RenameCells(object = Integration2, new.names = SCTIDs)
Integration2@meta.data$IntegrationCluster <- Idents(Integration2)
##Jaccard Index

head(Integration1[[]])
head(Integration2[[]])
#Jaccard Similarity Heatmap1
pdf(file=paste0(output,"_PairWiseJaccardSetsHeatmap_PubvNaive.pdf"), width=27, height=25)
PairWiseJaccardSetsHeatmap(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                           set_names(Integration1@meta.data$IntegrationCluster, nm=colnames(Integration1)),
                           show_row_dend = F, show_column_dend = F, cluster_row = F, cluster_column =F,
                           col_low = "#FFFFFF", col_high = "#004D40")
dev.off()

#Jaccard Index Calculation
PairWiseJaccardSets <- PairWiseJaccardSets(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                                           set_names(Integration1@meta.data$IntegrationCluster, nm=colnames(Integration1)))
write.csv(PairWiseJaccardSets, file = paste0(output, "_PairWiseJaccardSets_PubvNaive.csv"))

#Jaccard Similarity Heatmap2
pdf(file=paste0(output,"_PairWiseJaccardSetsHeatmap_PubvSCT.pdf"), width=27, height=25)
PairWiseJaccardSetsHeatmap(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                           set_names(Integration2@meta.data$IntegrationCluster, nm=colnames(Integration2)),
                           show_row_dend = F, show_column_dend = F, cluster_row = F, cluster_column =F,
                           col_low = "#FFFFFF", col_high = "#004D40")
dev.off()

#Jaccard Index Calculation
PairWiseJaccardSets <- PairWiseJaccardSets(set_names(Published@meta.data$PublishedCluster, nm=colnames(Published)),
                                           set_names(Integration2@meta.data$IntegrationCluster, nm=colnames(Integration2)))
write.csv(PairWiseJaccardSets, file = paste0(output, "_PairWiseJaccardSets_PubvSCT.csv"))

#Jaccard Similarity Heatmap3
pdf(file=paste0(output,"_PairWiseJaccardSetsHeatmap_NaivevSCT.pdf"), width=27, height=25)
PairWiseJaccardSetsHeatmap(set_names(Integration1@meta.data$IntegrationCluster, nm=colnames(Integration1)),
                           set_names(Integration2@meta.data$IntegrationCluster, nm=colnames(Integration2)),
                           show_row_dend = F, show_column_dend = F, cluster_row = F, cluster_column =F,
                           col_low = "#FFFFFF", col_high = "#004D40")
dev.off()

#Jaccard Index Calculation
PairWiseJaccardSets <- PairWiseJaccardSets(set_names(Integration1@meta.data$IntegrationCluster, nm=colnames(Integration1)),
                                           set_names(Integration2@meta.data$IntegrationCluster, nm=colnames(Integration2)))
write.csv(PairWiseJaccardSets, file = paste0(output, "_PairWiseJaccardSets_NaivevSCT.csv"))
