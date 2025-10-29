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
BNST <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered3.rds")


data <- Pct$data
data.df <- data.frame(data)
filt <- data.df %>% filter(pct.exp > 25)
genes <- filt %>% select(features.plot)
genes <- unique(genes)
GenesBNST <- unlist(genes)

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

Allen.data <- read.csv("Allen_Data/Subsampled_30K.csv", header=TRUE, row.names=1)
Allen.data <- t(Allen.data)
Allen <- CreateSeuratObject(counts=Allen.data, project="Allen_10x", min.cells=10, min.features=200)





U1 <- union(GenesBNST, GenesMeA)
U2 <- union(U1, GenesPOA)
UF <- union(U2, GenesVMH)
write.table(UF, file="loom/Enriched_TF_list.txt")

mySeurat <- merge(x=BNST, y=c(MeA,POA,VMH, Allen), add.cell.ids=c("BNST","MeA","POA","VMH","Allensub"))
mySeurat
#saveRDS(mySeurat, file="loom/Merged_Seurat_AllRegions_Preg.rds")
mySeurat.loom <- as.loom(mySeurat, filename="loom/All_Merged_Preg_With_30KOutgroup_corrected_filtered.loom")
mySeurat.loom$close_all()