#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POABNST_DEExclusionClusterTests/BNST_DEincluded/vln/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_POA_PrimedvUnprimedClusterTests/BNST_PrimedvUnprimed.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/BNST_PrimedvUnprimed.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.hormone <- paste(Idents(mySeurat), mySeurat$hormone, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.hormone"
for (g in genelist){
pdf(file=paste0(output,g,".pdf"))
for (i in 0:43){
try({
ident1 <- paste0(i,"_Primed")
ident2 <- paste0(i, "_Unprimed")
idents <- c(ident1,ident2)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=1, combine=TRUE)
plot(vlnvar)
})
}
graphics.off()
}

