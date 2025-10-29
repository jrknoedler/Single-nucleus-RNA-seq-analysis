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
memory.limit(64000)

allen.data <- read.csv("Allen_Data/matrix.csv", header=TRUE, row.names=1)
allen.metadata <- read.csv("Allen_Data/metadata.csv", row.names=1, header=TRUE)
allen <- CreateSeuratObject(counts=allen.data, meta.data=allen.metadata)
head(allen[[]])
allen <- SetIdent(M, value=allen@meta.data$cluster_label)
allenidents <- levels(Idents(allen))

Master_allen = data.frame()

    for (a in Master_allen){(
      try({
        sub <- subset(Seurat, idents=c(a))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- rowSums(counts)
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_allen = rbind(Master_allen,Pseudo.df)
      })
    )}
write.csv(Master_allen, file="Allen_Data/Allen_Pseudobulk.csv")

Master_GRN = data.frame()
#Load Seurat objects and get detected genes
BNST <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered3.rds")
BNSTno <- nlevels(Idents(BNST))
newBNST.ids <- c(1:BNSTno)
names(newBNST.ids) <- levels(BNST)
BNST <- RenameIdents(BNST,newBNST.ids)
MeA <- readRDS("X:/Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_SCT.rds")
MeAno <- nlevels(Idents(MeA))
newMeA.ids <- c(1:MeAno)
names(newMeA.ids) <- levels(MeA)
MeA<- RenameIdents(MeA,newMeA.ids)
POA <- readRDS("Z:/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")
POAno <- nlevels(Idents(POA))
newPOA.ids <- c(1:POAno)
names(newPOA.ids) <- levels(POA)
POA<- RenameIdents(POA,newMeA.ids)
VMH <- readRDS("X:/Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")
VMHno <- nlevels(Idents(VMH))
newVMH.ids <- c(1:VMHno)
names(newVMH.ids) <- levels(VMH)
VMH <- RenameIdents(VMH,newVMH.ids)
Regions <- c(BNST, MeA, POA, VMH)
Conds <- c("Intact","Primed","Unprimed","Pregnant")


Master_VMH =data.frame()
for (c in Conds){
  try({
    Seurat <- subset(x=VMH, subset=Hormone==c)
    Total <- nlevels(Idents(Seurat))
    for (i in 1:Total){(
      try({
        sub <- subset(Seurat, idents=c(i))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- rowSums(counts)
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_VMH = rbind(Master_VMH,Pseudo.df)
      })
    )}
    
  })
}

Master_BNST =data.frame()
for (c in Conds){
  try({
    Seurat <- subset(x=BNST, subset=Hormone==c)
    Total <- nlevels(Idents(Seurat))
    for (i in 1:Total){(
      try({
        sub <- subset(Seurat, idents=c(i))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- rowSums(counts)
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_BNST = rbind(Master_BNST,Pseudo.df)
      })
    )}
    
  })
}

Master_MeA =data.frame()
for (c in Conds){
  try({
    Seurat <- subset(x=MeA, subset=Hormone==c)
    Total <- nlevels(Idents(Seurat))
    for (i in 1:Total){(
      try({
        sub <- subset(Seurat, idents=c(i))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- rowSums(counts)
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_MeA = rbind(Master_MeA,Pseudo.df)
      })
    )}
    
  })
}

Master_POA =data.frame()
for (c in Conds){
  try({
    Seurat <- subset(x=POA, subset=Hormone==c)
    Total <- nlevels(Idents(Seurat))
    for (i in 1:Total){(
      try({
        sub <- subset(Seurat, idents=c(i))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- rowSums(counts)
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_POA = rbind(Master_POA,Pseudo.df)
      })
    )}
    
  })
}




dim(Master_GRN)

Master_GRN <- t(Master_GRN)
