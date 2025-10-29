#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Guiping_Collaboration/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_MergedCCA_Filtered2_res1.2.rds")


data <- mySeurat@assays[["RNA"]]@counts
data <- as.data.frame(data)
data <- t(data)
idents <- mySeurat@active.ident
ident <- as.data.frame(idents)
write.table(ident, file=paste0(output, "MeAclusters.tsv"))
hormone <- mySeurat$Hormone
hormone <- as.data.frame(hormone)
write.table(hormone, file=paste0(output, "MeAcondition.tsv"))
metadata <- merge(ident, hormone, by="row.names")
write.table(metadata, file=paste0(output,"MeAmetadata.tsv"))
write.table(data, file=paste0(output,"MeA_countmatrix.tsv"))

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")


data <- mySeurat@assays[["RNA"]]@counts
data <- as.data.frame(data)
data <- t(data)
idents <- mySeurat@active.ident
ident <- as.data.frame(idents)
write.table(ident, file=paste0(output, "POAclusters.tsv"))
hormone <- mySeurat$Hormone
hormone <- as.data.frame(hormone)
write.table(hormone, file=paste0(output, "POAcondition.tsv"))
write.table(metadata, file=paste0(output,"POAmetadata.tsv"))
write.table(data, file=paste0(output,"POA_countmatrix.tsv"))

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")


data <- mySeurat@assays[["RNA"]]@counts
data <- as.data.frame(data)
data <- t(data)
idents <- mySeurat@active.ident
ident <- as.data.frame(idents)
write.table(ident, file=paste0(output, "VMHclusters.tsv"))
hormone <- mySeurat$Hormone
hormone <- as.data.frame(hormone)
write.table(hormone, file=paste0(output, "VMHcondition.tsv"))
write.table(metadata, file=paste0(output,"VMHmetadata.tsv"))
write.table(data, file=paste0(output,"VMH_countmatrix.tsv"))

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered2_sexclude_malat1regress_30pcsres1.2FINAL.rds")


data <- mySeurat@assays[["RNA"]]@counts
data <- as.data.frame(data)
data <- t(data)
idents <- mySeurat@active.ident
ident <- as.data.frame(idents)
write.table(ident, file=paste0(output, "BNSTclusters.tsv"))
hormone <- mySeurat$Hormone
hormone <- as.data.frame(hormone)
write.table(hormone, file=paste0(output, "BNSTcondition.tsv"))
metadata <- merge(ident, hormone, by="row.names")
write.table(metadata, file=paste0(output,"BNSTmetadata.tsv"))
write.table(data, file=paste0(output,"BNST_countmatrix.tsv"))