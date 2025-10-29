#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_Kiss1exploratory"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40")
names(new.cluster.ids) <- levels(POA)
POA <- RenameIdents(POA, new.cluster.ids)


pdf(file=paste0(output,"Tacsignaling.pdf"), width=30, heigh=20)
VlnPlot(POA, features=c("Tac1","Tac2","Tacr1","Tacr3"), pt.size=0, ncol=1)
dev.off()


Kiss1 <- subset(POA, Kiss1 > 0)

Kiss1 <- SetIdent(Kiss1, value=Kiss1@meta.data$Hormone)
pdf(file=paste0(output, "Kiss1posonlyl"), height=15)
VlnPlot(Kiss1,  features=c("Kiss1","Tac1","Tac2","Tacr1","Tacr3"), pt.size=0, ncol=1)
dev.off()