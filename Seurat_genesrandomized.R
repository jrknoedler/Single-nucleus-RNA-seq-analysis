#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_Doubletsremoved/RandomGenes"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
DefaultAssay(mySeurat) <- "RNA"
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
random.genes <- sample(genes.10x, 344)
mySeurat$celltype.sex <- paste(Idents(mySeurat), mySeurat$sex, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.sex"
for (i in 0:38){
try({
ident1 <- paste0(i,"_Male")
ident2 <- paste0(i,"_Female")
sex.dimorphism <- FindMarkers(mySeurat, ident.1 = ident1, ident.2=ident2, min.pct=0, logfc.threshold=0,  verbose=TRUE, features=random.genes)
write.csv(sex.dimorphism, file=paste0(output,i,"_permutated1.csv"))
})
}
