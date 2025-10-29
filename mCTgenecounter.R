#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MetaDEGs/MvP_AllDEGs_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumlighterfiltered_sexclude_malat1regress_30pcsres1.2.rds")
coding <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Proteincoding_Genesybols.txt")
coding <- unlist(coding)
coding <- unique(coding)
head(coding)
length(coding)
dim(coding)
mySeurat <- subset(mySeurat, idents=c("0"))

mat <- mySeurat@assays[["RNA"]]@counts

df <- as.data.frame(mat)
dim(df)
collapsed <- rowSums(df)
head(collapsed)
dim(collapsed)
length(collapsed[collapsed>0])
genes <- collapsed[collapsed >0]
length(genes)
genes <- as.data.frame(genes)
dim(genes)
head(genes)
genes <- rownames(genes)
head(genes)
genes <- noquote(genes)
head(genes)
final <- intersect(genes,coding)
length(final)
