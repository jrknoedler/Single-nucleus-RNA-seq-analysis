#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/VMH_IndependentAnalysis/VMH_Merged_Clusterfuck_Marchredo_25pct_"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/VMH_MvP_1.5.txt")
eDEG <- read.table("Genelists/VMH_PvU_1.5.txt")
sDEG2 <- read.table("Genelists/VMH_MvU_1.5.txt")
unlist(eDEG)
unlist(sDEG2)
unlist(Dimorphic)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered <- as.matrix(eDEG.filtered)

sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered <- as.matrix(sDEG2.filtered)

Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)
mySeurat <- NormalizeData(mySeurat, normalization.method="RC", scale.factor=1e6)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)


Pct <- DotPlot(mySeurat, features=Dimorphic.filtered)
head(Pct)
Pct[[]]
data <- Pct$data
head(data)