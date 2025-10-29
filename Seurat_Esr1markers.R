#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimedMerged_Exploratory/POA_MalevUnprimed_Esr1Markers_sexsplit"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/POA_MalevUnprimed_Esr1filtered.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output, "_set1.pdf"))
VlnPlot(mySeurat, features=c("Npy1r","Lhx9","Fbn2"), pt.size=0, split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "set2.pdf"))
VlnPlot(mySeurat, features=c("Nrip1","Ecel1","Greb1"),pt.size=0, split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "set3.pdf"))
VlnPlot(mySeurat, features=c("Tshz2","Garem1","Moxd1"), pt.size=0, split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "set4.pdf"))
VlnPlot(mySeurat, features=c("Lhx6","Zfp536","Crim1"), pt.size=0, split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "set5.pdf"))
VlnPlot(mySeurat, features=c("Nfia","Nfib","Nfix"), pt.size=0, split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "set6.pdf"))
VlnPlot(mySeurat, features=c("Trhr","Pdzrn4","Lncenc1"), pt.size=0, split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "set7.pdf"))
VlnPlot(mySeurat, features=c("Rarb","Slc18a2","Lix1"), pt.size=0, split.by="sex", ncol=1)
dev.off()
pdf(paste0(output, "set8.pdf"))
VlnPlot(mySeurat, features=c("Sytl4","Olfm3","March11"), pt.size=0, split.by="sex", ncol=1)
dev.off()
