#!/usr/bin/env Rscript
library(Seurat)
library(rgeos)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)
library(SeuratDisk)
library(scclusteval)


Master_GRN = data.frame()
#Load Seurat objects and get detected genes
BNST <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered3.rds")

MeA <- readRDS("Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_SCT.rds")

POA <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")

VMH <- readRDS("Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")

Regions <- c(BNST, MeA, POA, VMH)
Conds <- c("Intact","Primed","Unprimed","Pregnant")

for (r in Regions){
try({
for (c in Conds){
try({
Seurat <- subset(x=r, idents="Hormone"==c)
Total <- nlevels(Idents(Seurat))
for (i in 0:Total){
try({
counts <- Seurat@assays[["RNA"]]@counts
counts <- t(counts)
Pseudobulk <- rowSums(counts)
Master_GRN = cbind(Pseudobulk, Master_GRN)
})
}
})
}
})
}
dim(Master_GRN)

Master_GRN <- t(Master_GRN)

GRN <- CreateSeuratObject(counts=Master_GRN, project="GRN, min.cells=0, min.features=0)

GRN

mySeurat.loom <- as.loom(mySeurat, filename="loom/Hypothalamus_Pseudobulk.loom")
mySeurat.loom$close_all()