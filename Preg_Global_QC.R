#!/usr/bin/env Rscript
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)

###load VMH
VMH <- readRDS("Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift_Filtered1.rds")
VMH$Region <- "VMH"

hormonestatusVMH <- VMH$Hormone
hormonestatusVMH.df <- data.frame(hormonestatusVMH)
hormonestatusVMH.df$Pub <- with(hormonestatusVMH.df, ifelse(hormonestatusVMH.df$hormonestatusVMH =="Pregnant", "NewPub", "OldPub"))
PubVMH <- hormonestatusVMH.df[,2]
VMH <- AddMetaData(VMH, PubVMH, col.name="Pub")

###load BNST
BNST <- readRDS("Seurat/BNST_Barcode_Transfer/BNST_Barcode_Lift_Filtered1.rds")
BNST$Region <- "BNST"

hormonestatusBNST <- BNST$Hormone
hormonestatusBNST.df <- data.frame(hormonestatusBNST)
hormonestatusBNST.df$Pub <- with(hormonestatusBNST.df, ifelse(hormonestatusBNST.df$hormonestatusBNST =="Pregnant", "NewPub", "OldPub"))
PubBNST <- hormonestatusBNST.df[,2]
BNST <- AddMetaData(BNST, PubBNST, col.name="Pub")

###load MeA
MeA <- readRDS("Seurat/MeA_Barcode_Transfer/MeA_Barcode_Lift_SCT.rds")
MeA$Region <- "MeA"

hormonestatusMeA <- MeA$Hormone
hormonestatusMeA.df <- data.frame(hormonestatusMeA)
hormonestatusMeA.df$Pub <- with(hormonestatusMeA.df, ifelse(hormonestatusMeA.df$hormonestatusMeA =="Pregnant", "NewPub", "OldPub"))
PubMeA <- hormonestatusMeA.df[,2]
MeA <- AddMetaData(MeA, PubMeA, col.name="Pub")

###Load POA
POA <- readRDS("/scratch/users/tsakura/Integration_Pregnancy/Naive_Integration_Fixed/POA_filtered_40/POA_Naive_Integration_Fixed_Filtered.rds")
POA$Region <- "POA"

hormonestatusPOA <- POA$Hormone
hormonestatusPOA.df <- data.frame(hormonestatusPOA)
hormonestatusPOA.df$Pub <- with(hormonestatusPOA.df, ifelse(hormonestatusPOA.df$hormonestatusPOA =="Pregnant", "NewPub", "OldPub"))
PubPOA <- hormonestatusPOA.df[,2]
POA <- AddMetaData(POA, PubPOA, col.name="Pub")

FullData <- merge(x=BNST, y=c(MeA,POA,VMH), add.cell.ids=c("BNST","MeA","POA","VMH"))
FullData[["percent.Malat1"]] <- PercentageFeatureSet(FullData, pattern = "Malat1")
FullData <- SCTransform(FullData, verbose=TRUE, vars.to.regress=c("orig.ident","percent.Malat1"))
FullData <- RunPCA(FullData)
FullData <- FindNeighbors(FullData, dims=1:50)
FullData <- FindClusters(FullData, resolution=1.5)
FullData <- RunUMAP(FullData, dims=1:50)
saveRDS(FullData, file="Pregnancy_Figures/All_snRNAseq_regressed.rds")
