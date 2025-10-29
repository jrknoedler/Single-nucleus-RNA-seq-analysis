#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_IndependentRound1_"


library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

POA <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")

POAdata <- POA@assays[["RNA"]]@counts
POAdata <- as.matrix(POAdata)
POAdata <- t(POAdata)
write.table(POAdata, file="geo_submission_8272021/BNSTcountmatrix.tsv",quote=FALSE, sep="\t", col.names=TRUE)
