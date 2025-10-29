#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_BusMergeTest/Vln_MvC/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/MeA_BusMergeTest/MeA_BusMergeTest_UMIcutoff_noUnprimed.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/MeA_MalevCastrate_3v2.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.status <- paste(Idents(mySeurat), mySeurat$status, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.status"
for (g in genelist){
pdf(file=paste0(output,g,".pdf"))
for (i in 0:46){
try({
ident1 <- paste0(i,"_IntactMale")
ident2 <- paste0(i,"_PrimedFemale")
ident3 <- paste0(i,"_CastrateMale")
idents <- c(ident1,ident2,ident3)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=1, combine=TRUE)
plot(vlnvar)
})
}
graphics.off()
}

