#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_Doubletsremoved/DiffgenesClustExp/UpMale/MalevUnprimed_UpinMaleTRAP_"

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
for (i in genelist){
try({
pdf(paste(output, i ,".pdf"))
Vln <- VlnPlot(mySeurat, features=c(i), pt.size=0)
print(Vln)
dev.off()
})
}

