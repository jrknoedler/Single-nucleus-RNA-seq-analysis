#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat//BNST_MalevUnprimedDiff/Aromatase/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_Aromataseneurons.rds")
DefaultAssay(mySeurat) <- "RNA"
head(mySeurat[[]])
genelist <- read.table("Genelists/BNST_MalevUnprimed.txt", header=FALSE)
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
clusters <- c("1","6","7","12","14","15","17","20")
for (i in clusters){
try({
ident1 <- paste0(i,"_Male")
ident2 <- paste0(i,"_Female")
sex.dimorphism <- FindMarkers(mySeurat, ident.1 = ident1, ident.2=ident2, min.pct=0, logfc.threshold=0,  verbose=TRUE, features=filtered.genelist, test.use="negbinom")
write.csv(sex.dimorphism, file=paste0(output,i,"_dimorphism.csv"))
})
}
