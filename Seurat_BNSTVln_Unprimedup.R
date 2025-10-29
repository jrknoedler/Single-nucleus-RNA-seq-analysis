#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_MalevUnprimedDiff/VlnUnprimedUP/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_MalevUnprimed_Exploratory/BNST_MaleUnprimed_doubletsremoved_filtered1.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/BNST_MalevUnprimed_UpUnprimed.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.sex <- paste(Idents(mySeurat), mySeurat$sex, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.sex"
for (g in genelist){
try(
dir.create(paste0(output,g)) 
)
for (i in 0:42){
try({
out.dir <- (paste0(output,"/",g,"/"))
pdf(paste0(out.dir,i,".pdf"))
ident1 <- paste0(i,"_Male")
ident2 <- paste0(i,"_Female")
idents <- c(ident1,ident2)
vln <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=1, combine=TRUE)
print(vln)
graphics.off()
}
)
}}

