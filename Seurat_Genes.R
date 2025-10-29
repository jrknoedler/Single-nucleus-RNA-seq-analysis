#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/Expressed_Genes/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)


PrimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentRound1__Primed.rds")
UnprimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentRound1__Unprimed.rds")
IntactPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentRound1__Intact.rds")
IntactMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentRound1__Intact.rds")

All.genes.PrimedVMH <- c(rownames(PrimedVMH))
write.csv(All.genes.PrimedVMH, file=paste0(output,"PrimedVMHgenes.csv"))

All.genes.UnprimedBNST <- c(rownames(UnprimedBNST))
write.csv(All.genes.UnprimedBNST, file=paste0(output,"UnprimedBNSTgenes.csv"))

All.genes.IntactPOA <- c(rownames(IntactPOA))
write.csv(All.genes.IntactPOA, file=paste0(output,"IntactPOAgenes.csv"))

All.genes.IntactMeA <- c(rownames(IntactMeA))
write.csv(All.genes.IntactMeA, file=paste0(output,"IntactMeAgenes.csv"))

