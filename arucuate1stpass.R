#!/usr/bin/env Rscript



output <- "Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_0.1CPM_hicat"
library(Seurat)
library(cowplot)
library(Matrix)
library(patchwork)
library(dplyr)
library(tidyverse)


mySeurat <- readRDS("arcuateSeurat_allmeta.rds")


mySeurat <- NormalizeData(mySeurat)
pdf(file="Arcuate_NRs.pdf", width=40, height=20)
VlnPlot(mySeurat, features=c("Esr1","Pgr"), ncol=1, pt.size=0, group.by="X7.clust_all")
dev.off()
pdf(file="Arcuate_NRs_dots.pdf", width=40, height=20)
VlnPlot(mySeurat, features=c("Esr1","Pgr"), ncol=1, pt.size=1, group.by="X7.clust_all")
dev.off()
pdf(file="Arcuate_NRs_Neurontypes.pdf", width=40, heigh=20)
VlnPlot(mySeurat, features=c("Esr1","Pgr"), ncol=1, pt.size=1, group.by="X8.clust_all_neurons")
dev.off()

pdf(file="Arcuate_NRs_Neurontypes2.pdf", width=40, heigh=20)
VlnPlot(mySeurat, features=c("Esr1","Pgr"), ncol=1, pt.size=1, group.by="X10.clust_neurons")
dev.off()



pdf(file="Arcuate_NRs_NeurontypesUMI.pdf", width=40, heigh=10)
VlnPlot(mySeurat, features=c("nCount_RNA"), ncol=1, pt.size=1, group.by="X10.clust_neurons")
dev.off()