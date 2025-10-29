#!/usr/bin/env Rscript

output <- "Seurat/POA_IndependentAnalysis/POA_Regulons"

library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(loomR)

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
pyScenicLoomFile <- file.path("pySCENIC_Singularity/POA_DefaultNES_Final/auc_mtx_filtered.loom")
loom <- open_loom(pyScenicLoomFile,mode="r")
exprMat <- get_dgem(loom)
cellInfo <- get_cellAnnotation(loom)
head(cellInfo)
regulons_incidMat <- get_regulons(loom, attrName='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulonsAuc(loom, attrName='RegulonsAUC')
regulonsAUC
regulonsAucThresholds <- get_regulonThresholds(loom)
regulonsAucThresholds
embeddings <- get_embeddings(loom)
embeddings
close_loom(loom)

regulonsAUC <- regulonsAUC[onlyNonDuplicatedExtended(rownames(regulonsAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(regulonsAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"allregulonheatmap.pdf"))
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()
Degulators <- read.table("RegulonsCompiled/POA/Degulators.txt")
Degulators <- unlist(Degulators)
Degulators <- as.matrix(Degulators)
Degulators <- regulonsAUC[Degulators]

Degulons <- read.table("RegulonsCompiled/POA/POA_SCENICDegulons.txt")
Degulons <- unlist(Degulons)
Degulons <- as.matrix(Degulons)
Degulons <- regulonsAUC[Degulons]

DegulatorActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(Degulators)[,cells]))
DegulatorActivity_byCellType_Scaled <- t(scale(t(DegulatorActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"Degulatorheatmap.pdf"))
pheatmap::pheatmap(DegulatorActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

DegulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(Degulons)[,cells]))
DegulonActivity_byCellType_Scaled <- t(scale(t(DegulonActivity_byCellType), center = T, scale=T))

pdf(file=paste0(output,"Degulonheatmap.pdf"))
pheatmap::pheatmap(DegulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=viridis(2000),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

