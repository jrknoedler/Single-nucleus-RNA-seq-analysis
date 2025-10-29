#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/PaperQC/PaperQC"
regoutput <- "RegulonsCompiled/POA/Sexmarkersperclust/"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")
MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")

MalePOA <- subset(x=POA, subset=Hormone=="Intact")
PrimedPOA <- subset(x=POA, subset=Hormone=="Primed")
UnprimedPOA <- subset(x=POA, subset=Hormone=="Unprimed")

MalePOA$region <- "POA"
MalePOA$sample <- "IntactPOA"
PrimedPOA$region <- "POA"
PrimedPOA$sample <- "PrimedPOA"
UnprimedPOA$region <- "POA"
UnprimedPOA$sample <- "UnprimedPOA"
MalePOA$batch <- "1"
PrimedPOA$batch <- "2"
UnprimedPOA$batch <- "3"

MaleBNST <- subset(x=BNST, subset=Hormone=="Intact")
PrimedBNST <- subset(x=BNST, subset=Hormone=="Primed")
UnprimedBNST <- subset(x=BNST, subset=Hormone=="Unprimed")

MaleBNST$region <- "BNST"
MaleBNST$sample <- "MaleBNST"
PrimedBNST$region <- "BNST"
PrimedBNST$sample <- "PrimedBNST"
UnprimedBNST$region <- "BNST"
UnprimedBNST$sample <- "UnprimedBNST"
MaleBNST$batch <- "4"
PrimedBNST$batch <- "2"
UnprimedBNST$batch <- "3"

MaleVMH <- subset(x=VMH, subset=Hormone=="Intact")
PrimedVMH <- subset(x=VMH, subset=Hormone=="Primed")
UnprimedVMH <- subset(x=VMH, subset=Hormone=="Unprimed")

MaleVMH$region <- "VMH"
MaleVMH$sample <- "MaleVMH"
PrimedVMH$region <- "VMH"
PrimedVMH$sample <- "PrimedVMH"
UnprimedVMH$region <- "VMH"
UnprimedVMH$sample <- "UnprimedVMH"
MaleVMH$batch <- "5"
PrimedVMH$batch <- "6"
UnprimedVMH$batch <- "7"


MaleMeA <- subset(x=MeA, subset=Hormone=="Intact")
PrimedMeA <- subset(x=MeA, subset=Hormone=="Primed")
UnprimedMeA <- subset(x=MeA, subset=Hormone=="Unprimed")

MaleMeA$region <- "MeA"
MaleMeA$sample <- "MaleMeA"
PrimedMeA$region <- "MeA"
PrimedMeA$sample <- "PrimedMeA"
UnprimedMeA$region <- "MeA"
UnprimedMeA$sample <- "UnprimedMeA"
MaleMeA$batch <- "5"
PrimedMeA$batch <- "6"
UnprimedMeA$batch <- "7"

mySeurat <- merge(MalePOA, y=c(PrimedPOA,UnprimedPOA,MaleBNST,PrimedBNST,UnprimedBNST,MaleVMH,PrimedVMH,UnprimedVMH,MaleMeA,PrimedMeA,UnprimedMeA), project="QC")
mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("sample","percent.Malat1"))
saveRDS(mySeurat, file=paste0(output, "prePCA.rds"))
pdf(file=paste0(output,"QCUMIplot.pdf"),width=15)
VlnPlot(mySeurat, features=c("nCount_RNA"), split.by="sample")
dev.off()
pdf(file=paste0(output,"QCgeneplot.pdf"),width=15)
VlnPlot(mySeurat, features=c("nFeature_RNA"), split.by="sample")
dev.off()
mySeurat <- RunPCA(mySeurat)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dims=1:50)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:50)
mySeurat <- FindClusters(mySeurat, resolution=2)
pdf(paste0(output,"_UMAP.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds"))
pdf(paste0(output,"_UMAP_sexlabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sex")
dev.off()

pdf(paste0(output,"_UMAP_hormonelabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()


pdf(paste0(output,"_UMAP_regionlabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="region")
dev.off()

pdf(paste0(output,"_UMAP_batchlabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="batch")
dev.off()


pdf(paste0(output,"_UMAP_idabel.pdf"), width=15,height=15)
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sample")
dev.off()