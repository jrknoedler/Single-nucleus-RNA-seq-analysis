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

Master_GRN = data.frame()
#Load Seurat objects and get detected genes
BNST <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered3.rds")
BNSTno <- nlevels(Idents(BNST))
newBNST.ids <- c(1:BNSTno)
names(newBNST.ids) <- levels(BNST)
BNST <- RenameIdents(BNST,newBNST.ids)
MeA <- readRDS("Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_SCT.rds")
MeAno <- nlevels(Idents(MeA))
newMeA.ids <- c(1:MeAno)
names(newMeA.ids) <- levels(MeA)
MeA<- RenameIdents(MeA,newMeA.ids)
POA <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")
POAno <- nlevels(Idents(POA))
newPOA.ids <- c(1:POAno)
names(newPOA.ids) <- levels(POA)
POA<- RenameIdents(POA,newMeA.ids)
VMH <- readRDS("Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")
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
        Pseudobulk <- apply(counts, 1, function(x)mean(x, trim=0.25))
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_VMH = rbind(Master_VMH,Pseudo.df)
      })
    )}
    
  })
}
remove(VMH)
write.csv(Master_VMH, file="loom/VMH_GRN_trimmed_pseudobulk.csv")
Master_BNST =data.frame()
for (c in Conds){
  try({
    Seurat <- subset(x=BNST, subset=Hormone==c)
    Total <- nlevels(Idents(Seurat))
    for (i in 1:Total){(
      try({
        sub <- subset(Seurat, idents=c(i))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- apply(counts, 1, function(x)mean(x, trim=0.25))
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_BNST = rbind(Master_BNST,Pseudo.df)
      })
    )}
    
  })
}
remove(BNST)
write.csv(Master_BNST, file="loom/BNST_GRN_trimmed_pseudobulk.csv")
Master_MeA =data.frame()
for (c in Conds){
  try({
    Seurat <- subset(x=MeA, subset=Hormone==c)
    Total <- nlevels(Idents(Seurat))
    for (i in 1:Total){(
      try({
        sub <- subset(Seurat, idents=c(i))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- apply(counts, 1, function(x)mean(x, trim=0.25))
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_MeA = rbind(Master_MeA,Pseudo.df)
      })
    )}
    
  })
}
remove(MeA)
write.csv(Master_MeA, file="loom/MeA_GRN_trimmed_pseudobulk.csv")
Master_POA =data.frame()
for (c in Conds){
  try({
    Seurat <- subset(x=POA, subset=Hormone==c)
    Total <- nlevels(Idents(Seurat))
    for (i in 1:Total){(
      try({
        sub <- subset(Seurat, idents=c(i))
        counts <- sub@assays[["RNA"]]@counts
        Pseudobulk <- apply(counts, 1, function(x)mean(x, trim=0.25))
        Pseudo.df <- data.frame(Pseudobulk)
        Pseudo.df <- t(Pseudo.df)
        Master_POA = rbind(Master_POA,Pseudo.df)
      })
    )}
    
  })
}

remove(POA)
write.csv(Master_POA, file="loom/POA_GRN_trimmed_pseudobulk.csv")
