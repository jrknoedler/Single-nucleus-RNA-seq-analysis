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

TFs <- read.table("/home/groups/nirao/pySCENIC/resources/mm_mgi_tfs.txt", header=FALSE)
TFs <- unlist(TFs)
#Load Seurat objects and get detected genes
#BNST <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered3.rds")
#DefaultAssay(BNST) <- "RNA"
#BNST$Region <- "BNST"
#genes.BNST <- (x=rownames(x=BNST))
#unlist(genes.BNST)
#expTB <- intersect(genes.BNST, TFs)
#Pct <- DotPlot(BNST, features=expTB)

#data <- Pct$data
#data.df <- data.frame(data)
#filt <- data.df %>% filter(pct.exp > 25)
#genes <- filt %>% select(features.plot)
#genes <- unique(genes)
#GenesBNST <- unlist(genes)
#write.csv(GenesBNST, file="loom/BNST_Tfs.txt")

MeA <- readRDS("Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_SCT.rds")
DefaultAssay(MeA) <- "RNA"
MeA$Region <- "MeA"
genes.MeA <- (x=rownames(x=MeA))
unlist(genes.MeA)
expTM <- intersect(genes.MeA, TFs)
Pct <- DotPlot(MeA, features=expTM)

data <- Pct$data
data.df <- data.frame(data)
filt <- data.df %>% filter(pct.exp > 25)
genes <- filt %>% select(features.plot)
genes <- unique(genes)
GenesMeA <- unlist(genes)
write.csv(GenesMeA, file="loom/MeA_Tfs.txt")

POA <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")
DefaultAssay(POA) <- "RNA"
POA$Region <- "POA"
genes.POA <- (x=rownames(x=POA))
unlist(genes.POA)
expTP <- intersect(genes.POA, TFs)
Pct <- DotPlot(POA, features=expTP)

data <- Pct$data
data.df <- data.frame(data)
filt <- data.df %>% filter(pct.exp > 25)
genes <- filt %>% select(features.plot)
genes <- unique(genes)
GenesPOA <- unlist(genes)
write.csv(GenesPOA, file="loom/POA_Tfs.txt")

VMH <- readRDS("Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")
DefaultAssay(VMH) <- "RNA"
VMH$Region <- "VMH"
genes.VMH <- (x=rownames(x=VMH))
unlist(genes.VMH)
expTV <- intersect(genes.VMH, TFs)
Pct <- DotPlot(VMH, features=expTV)

data <- Pct$data
data.df <- data.frame(data)
filt <- data.df %>% filter(pct.exp > 25)
genes <- filt %>% select(features.plot)
genes <- unique(genes)
GenesVMH <- unlist(genes)
write.csv(GenesVMH, file="loom/VMH_Tfs.txt")

Allen.data <- read.csv("Allen_Data/Subsampled_30K.csv", header=TRUE, row.names=1)
Allen.data <- t(Allen.data)
Allen <- CreateSeuratObject(counts=Allen.data, project="Allen_10x", min.cells=50, min.features=200)


#BNSTAllen <- merge(x=BNST, y=Allen, add.cell.ids=c("BNST","Allensub"))
#saveRDS(BNSTAllen, file="loom/Merged_BNST_Allen_Preg.rds")
#BNSTAllen.loom <- as.loom(BNSTAllen, filename="loom/BNST_Allen.loom")
#BNSTAllen.loom$close_all()
MeAAllen <- merge(x=MeA, y=Allen, add.cell.ids=c("MeA","Allensub"))
saveRDS(MeAAllen, file="loom/Merged_MeA_Allen_Preg.rds")
MeAAllen.loom <- as.loom(MeAAllen, filename="loom/MeA_Allen.loom")
MeAAllen.loom$close_all()
POAAllen <- merge(x=POA, y=Allen, add.cell.ids=c("POA","Allensub"))
saveRDS(POAAllen, file="loom/Merged_POA_Allen_Preg.rds")
POAAllen.loom <- as.loom(POAAllen, filename="loom/POA_Allen.loom")
POAAllen.loom$close_all()
VMHAllen <- merge(x=VMH, y=Allen, add.cell.ids=c("VMH","Allensub"))
saveRDS(VMHAllen, file="loom/Merged_VMH_Allen_Preg.rds")
VMHAllen.loom <- as.loom(VMHAllen, filename="loom/VMH_Allen.loom")
VMHAllen.loom$close_all()

