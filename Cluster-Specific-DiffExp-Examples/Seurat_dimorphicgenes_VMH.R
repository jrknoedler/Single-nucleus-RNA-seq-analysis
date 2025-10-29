#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_MalevPrimed_SCTDE/StandardTest/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(scater)
mySeurat <- readRDS("Seurat/VMH_MalevPrimedMerged_Exploratory/VMH_MalevPrimed_Exploratory_filtered2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
head(mySeurat[[]])
genelist <- read.table("Genelists/VMH_MalevPrimedAll.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
mySeurat$celltype.sex <- paste(Idents(mySeurat), mySeurat$sex, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.sex"
for (i in 0:33){
try({
ident1 <- paste0(i,"_Male")
ident2 <- paste0(i,"_Female")
sex.dimorphism <- FindMarkers(mySeurat, ident.1 = ident1, ident.2=ident2, min.pct=0, logfc.threshold=0.1,  verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(output,i,"_allgenes_dimorphism.csv"))
})
}
