#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Vln_Fig6/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("RNASeqKmeans/RNASeqKmeans/BNST_1.5unique.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.Hormone <- paste(Idents(mySeurat), mySeurat$Hormone, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.Hormone"
for (g in genelist){
pdf(file=paste0(output,g,".pdf"))
for (i in 0:33){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i, "_Primed")
ident3 <- paste0(i, "_Unprimed")
idents <- c(ident1,ident2,ident3)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=0, combine=TRUE, cols=c("blue", "red", "green"))
plot(vlnvar)
})
}
graphics.off()
}

