#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_CelltypesDraft_Malat1regress_Esr1filtered_sexclude_res1.2alt2reordered"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(file=paste0(output,"Tacr1_peptides.pdf"), height=20, width=25)
VlnPlot(mySeurat, features=c("Tacr1","Tac1","Cartpt","Cck","Sst","Vip","Tac2","Tac3"), ncol=1, pt.size=0)
dev.off()
