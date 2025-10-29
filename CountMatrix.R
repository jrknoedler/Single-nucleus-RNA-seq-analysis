#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_IndependentRound1_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")

POAdata <- POA@assays[["RNA"]]@counts
POAdata <- as.matrix(POAdata)
POAdata <- t(POAdata)
write.table(POAdata, file="geo_submission_8272021/POAcountmatrix.tsv",quote=FALSE, sep="\t", col.names=TRUE)


MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")

MeAdata <- MeA@assays[["RNA"]]@counts
MeAdata <- as.matrix(MeAdata)
MeAdata <- t(MeAdata)
write.table(MeAdata, file="geo_submission_8272021/MeAcountmatrix.tsv",quote=FALSE, sep="\t", col.names=TRUE)


VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
VMHdata <- VMH@assays[["RNA"]]@counts
VMHdata <- as.matrix(VMHdata)
VMHdata <- t(VMHdata)
write.table(VMHdata, file="geo_submission_8272021/VMHcountmatrix.tsv",quote=FALSE, sep="\t", col.names=TRUE)

