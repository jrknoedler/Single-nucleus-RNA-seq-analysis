#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_Doubletsremoved/DimorphicTest/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(scater)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
head(mySeurat[[]])
genelist <- read.table("Genelists/POA_MalevUnprimed.txt", header=FALSE)
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
for (i in 0:38){
try({
ident1 <- paste0(i,"_Male")
ident2 <- paste0(i,"_Female")
sex.dimorphism <- FindMarkers(mySeurat, ident.1 = ident1, ident.2=ident2, min.pct=0, logfc.threshold=0.1,  verbose=TRUE, test.use="MAST")
write.csv(sex.dimorphism, file=paste0(output,i,"_MAST_allgenes_dimorphism.csv"))
})
}
