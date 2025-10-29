#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_MalevPrimed_StandardCCA/Vln/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/VMH_MalevPrimed_StandardCCA/VMH_MalevPrimed_StandardCCA_.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/VMH_MalevPrimedAll.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.sex <- paste(Idents(mySeurat), mySeurat$sex, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.sex"
for (g in genelist){
pdf(file=paste0(output,g,".pdf"))
for (i in 0:29){
try({
ident1 <- paste0(i,"_Male")
ident2 <- paste0(i,"_Female")
idents <- c(ident1,ident2)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=1, combine=TRUE)
plot(vlnvar)
})
}
graphics.off()
}

