#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_BusMergeTest/UnfilteredVln/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_BusMergeTest/POA_BusMergeTest.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/POA_Allcomps_new.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.status <- paste(Idents(mySeurat), mySeurat$status, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.status"
for (g in genelist){
pdf(file=paste0(output,g,".pdf"))
for (i in 0:52){
try({
ident1 <- paste0(i,"_IntactMale")
ident2 <- paste0(i,"_PrimedFemale")
ident3 <- paste0(i,"_UnprimedFemale")
ident4 <- paste0(i,"_CastrateMale")
idents <- c(ident1,ident2,ident3,ident4)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=1, combine=TRUE)
plot(vlnvar)
})
}
graphics.off()
}

