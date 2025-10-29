#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/Paper1stpass_"
library(liger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(cowplot)
library(Matrix)
library(patchwork)
library(dplyr)
library(tidyverse)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)

datar <- round(data)

RoundedSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone

RoundedSeurat <- AddMetaData(mySeurat, metadata=idents, col.name="Ident")
RoundedSeurat <- AddMetaData(mySeurat, metadata=hormone, col.name="Hormone")
Idents(RoundedSeurat) <- "Ident"
head(RoundedSeurat[[]])
