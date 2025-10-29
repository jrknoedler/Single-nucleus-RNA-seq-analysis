#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimed_Doubletsremoved/POA_MalevUnprimed_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
mySeurat
mySeurat[[]]
head(mySeurat[[]])
antisense <- read.table("Genelists/antisensenames.txt", header=FALSE)
antisense <- as.matrix(antisense)
head(antisense)
mySeurat[["percent.antisense"]] <- PercentageFeatureSet(mySeurat, features=antisense)
pdf(file=paste0(output, "_antisensecounts.pdf"))
VlnPlot(mySeurat, features = c("percent.antisense"), split.by="sex"))
dev.off()
