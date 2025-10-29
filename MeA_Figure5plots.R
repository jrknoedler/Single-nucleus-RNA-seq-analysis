#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/MeA_ngf_Fig5"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)



mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")
SFARI <- read.table("

pdf(file=paste0(output,"UMAP_Hormonelabelshuffle.pdf"))
DimPlot(mySeurat, reduction="umap", group.by="Hormone", shuffle=TRUE, seed=1, label=FALSE, pt.size=0.4, cols=c("blue", "red", "green"))
dev.off()
pdf(file=paste0(output,"brainstormvln.pdf"), width=40, height=80)
VlnPlot(mySeurat, features=c("Esr1","Ar","Gad1","Slc17a6","Slc18a2","Sst","Cartpt","Gal","Cck","Esr2","Tac1","Tac2","Th","Maob","Kiss1","Pappa","Col25a1","Npy","Fbn2"), ncol=1, pt.size=0)
dev.off()
ngf <- subset(mySeurat, idents=c("4","9","25"))
ngf.markers <- FindAllMarkers(ngf, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.50, test.use="MAST", min.diff.pct = 0)
write.csv(ngf.markers, file=paste0(output,"ngfmarkers.csv"))