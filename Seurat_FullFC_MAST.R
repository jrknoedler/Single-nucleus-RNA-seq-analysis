#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_BusMergeTest/BNSTFullFCtest_renamed_MAST"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/BNST_BusMergeTest/BNST_BusMergeTest.rds")
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
mySeurat.markers1 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="IntactMale", ident.2="PrimedFemale", min.pct=0, logfc.threshold = 0, test.use="MAST")
write.csv(mySeurat.markers1, file=paste0(output,"_MalevPrimed.csv"))
mySeurat.markers2 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="PrimedFemale", ident.2="IntactMale", min.pct=0, logfc.threshold = 0, test.use="MAST")
write.csv(mySeurat.markers2, file=paste0(output,"_PrimedvMale.csv"))
mySeurat.markers3 <- FindMarkers(mySeurat, features=filtered.genelist, ident.1="PrimedFemale", ident.2="UnprimedFemale", min.pct=0, logfc.threshold = 0, test.use="MAST")
write.csv(mySeurat.markers3, file=paste0(output,"_PrimedvUnprimed.csv"))
mySeurat.markers4 <- FindAllMarkers(mySeurat, features=filtered.genelist, ident.1="IntactMale", ident.2="UnprimedFemale", min.pct=0, logfc.threshold = 0, test.use="MAST")
write.csv(mySeurat.markers4, file=paste0(output,"_MalevUnprimed.csv"))
