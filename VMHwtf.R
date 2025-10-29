#!/usr/bin/env Rscript

output <- "Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_20pcs"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
PrimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1_Primed_30pcs_res1_Primed.rds")
MaleVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1_Intact_30pcs_res1_Primed.rds")
UnprimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered3_Unprimed_30pcs_res1.rds")
#PrimedVMH <- SCTransform(PrimedVMH)
#MaleVMH <- SCTransform(MaleVMH)
#UnprimedVMH <- SCTransform(UnprimedVMH)
MaleVMH$sex <- "Male"
MaleVMH$Hormone <- "Intact"
PrimedVMH$sex <- "Female"
PrimedVMH$Hormone <- "Primed"
UnprimedVMH$sex <- "Female"
UnprimedVMH$Hormone <- "Unprimed"
mySeurat <- merge(MaleVMH, y=c(PrimedVMH,UnprimedVMH), add.cell.ids=c("MaleVMH","PrimedVMH","UnprimedVMH"), project="Merged")
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress="orig.ident")
mySeurat <- RunPCA(mySeurat)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:20)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:20)
mySeurat <- FindClusters(mySeurat, resolution=1)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))