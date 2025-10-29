#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_BusMergeTest/BNSTFullFCtest_mvp_cDNA"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/BNST_BusMergeTest/BNST_BusMergeTest_MvP_Hailmary_cDNAonly.rds")
mySeurat <- subset(mySeurat, idents=c("8"), invert=TRUE)
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
mySeurat$celltype.status <- paste(mySeurat$status)
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.status"
Idents(mySeurat)
genelist <- read.table("Genelists/BNST_allgenes.txt", header=FALSE)
genelist <- unlist(genelist)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
mySeurat.markers1 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="IntactMale", ident.2="PrimedFemale", min.pct=0, logfc.threshold = 0)
write.csv(mySeurat.markers1, file=paste0(output,"_MalevPrimed.csv"))

