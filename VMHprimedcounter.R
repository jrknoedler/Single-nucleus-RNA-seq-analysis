#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_Merge_Fig6storyboard_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)
mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")

subset1 <- subset(x=mySeurat, subset=Hormone=="Primed")

prop.table(table(Idents(subset1))) 
head(subset1[[]])
excit <- subset(mySeurat, idents=c("1","2","3","4","6","8","14","16","17","19","20","21","23","24","25","26"))
excitprimed <- subset(x=excit, subset=Hormone=="Primed")
prop.table(table(Idents(excitprimed)))

mySeurat2 <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered3_Primed_10pcs_res1_Primed.rds")
prop.table(table(Idents(mySeurat2)))

Filt <- subset(mySeurat, idents=c("9","10"), invert=TRUE)
Filt <- subset(Filt, subset=Hormone=="Primed")
prop.table(table(Idents(Filt)))


vl <- subset(mySeurat, idents=c("6","20","23","25"), invert=TRUE)
vl <- subset(vl, subset=Hormone=="Primed")
prop.table(table(Idents(vl)))

CCA <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_CCA20pcs.rds")
prop.table(table(Idents(CCA), CCA$Hormone), margin=2)

pdf(file=paste0(output,"CCA_Cckar_Hormonesplit.pdf"), width=30)
VlnPlot(CCA, features=c("Cckar"), pt.size=0, split.by="Hormone")
dev.off()
pdf(file=paste0(output,"CCA_Esr1.pdf"), width=30)
VlnPlot(CCA, features=c("Esr1"), pt.size=0)
dev.off()
ccafilt <- readRDS("Seurat/VMH_IndependentAnalysis/VMHCCAfiltered1.rds")
prop.table(table(Idents(ccafilt), ccafilt$Hormone), margin=2)
pdf(file=paste0(output,"CCA_Cckafiltr_Hormonesplit.pdf"), width=30)
VlnPlot(ccafilt, features=c("Cckar"), pt.size=0, split.by="Hormone")
dev.off()
pdf(file=paste0(output,"CCA_filtr_Vglut.pdf"), width=30)
VlnPlot(ccafilt, features=c("Slc17a6"), pt.size=0)
dev.off()
pdf(file=paste0(output,"CCA_Esr1filtr_Hormonesplit.pdf"), width=30)
VlnPlot(ccafilt, features=c("Esr1"), pt.size=0, split.by="Hormone")
dev.off()
filtexcit <- subset(ccafilt, idents=c("0","1","4","6","7","13","23","24","28","30","31"), invert=TRUE)
prop.table(table(Idents(filtexcit), filtexcit$Hormone), margin=2)