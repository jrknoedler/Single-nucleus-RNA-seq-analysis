#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_BusMergeTest/POAFullFCtest2_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/POA_BusMergeTest/POA_BusMergeTest_filtered1.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
mySeurat$celltype.status <- paste(mySeurat$status)
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.status"
Idents(mySeurat)
genelist <- read.table("Genelists/VMH_allgenes.txt", header=FALSE)
genelist <- unlist(genelist)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
mySeurat.markers1 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="IntactMale", ident.2="PrimedFemale", min.pct=0, logfc.threshold = 0)
write.csv(mySeurat.markers1, file=paste0(output,"_MalevPrimed.csv"))
mySeurat.markers2 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="PrimedFemale", ident.2="UnprimedFemale", min.pct=0, logfc.threshold = 0)
write.csv(mySeurat.markers2, file=paste0(output,"_PrimedvUnprimed.csv"))
mySeurat.markers4 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="IntactMale", ident.2="UnprimedFemale", min.pct=0, logfc.threshold = 0)
write.csv(mySeurat.markers4, file=paste0(output,"_MalevUnprimed.csv"))
mySeurat.markers3 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="IntactMale", ident.2="CastrateMale", min.pct=0, logfc.threshold = 0)
write.csv(mySeurat.markers3, file=paste0(output,"_MalevCastrate.csv"))
