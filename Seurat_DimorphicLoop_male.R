#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_Doubletsremoved/DiffgenesClustExp/UpMale/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/POA_UpMalevUnprimed.txt", header=FALSE)
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
for (i in 0:38){
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

