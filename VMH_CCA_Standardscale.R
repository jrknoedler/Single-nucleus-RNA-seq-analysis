output <- "Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_CCAtestallstandardscale"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)


PrimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1_Primed_30pcs_res1_Primed.rds")
MaleVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1_Intact_30pcs_res1_Primed.rds")
UnprimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered3_Unprimed_30pcs_res1.rds")
list <- c(MaleVMH,UnprimedVMH, PrimedVMH )
genes.male <- (x=rownames(x=MaleVMH))
genes.unprimed <- (x=rownames(x=UnprimedVMH))
genes.primed <- (x=rownames(x=PrimedVMH))

allgenes <- intersect(intersect(genes.male,genes.primed),genes.unprimed)
DEGs <- read.csv("topGO/Total_ByRegion/VMH_Genesonly.txt", header=FALSE)
DEGs <- unlist(DEGs)
DEGs <- unique(DEGs)
DEGs.filtered <- intersect(DEGs, allgenes)
genelist <- read.table("ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)

DefaultAssay(PrimedVMH) <- "RNA"
DefaultAssay(UnprimedVMH) <- "RNA"
DefaultAssay(MaleVMH) <- "RNA"

PrimedVMH <- NormalizeData(PrimedVMH)
UnprimedVMH <- NormalizeData(UnprimedVMH)
MaleVMH <- NormalizeData(MaleVMH)

PrimedVMH <- FindVariableFeatures(PrimedVMH, selection.method="vst", nfeatures=2000)

UnprimedVMH <- FindVariableFeatures(UnprimedVMH, selection.method="vst", nfeatures=2000)

MaleVMH <- FindVariableFeatures(MaleVMH, selection.method="vst", nfeatures=2000)

features <- SelectIntegrationFeatures(object.list=list)

features.final <- union(features, DEGs.filtered)
features.final < unlist(features.final)
VMH.anchors <- FindIntegrationAnchors(object.list=list, anchor.features = features.final)

VMH.comb <- IntegrateData(anchorset=VMH.anchors)

DefaultAssay(VMH.comb) <- "integrated"

VMH.comb <- ScaleData(VMH.comb)
VMH.comb <- RunPCA(VMH.comb)
VMH.comb <- RunUMAP(VMH.comb, reduction="pca", dims=1:30)
VMH.comb <- FindNeighbors(VMH.comb, reduction="pca", dims=1:30)
VMH.comb <- FindClusters(VMH.comb, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(VMH.comb, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(VMH.comb, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(VMH.comb, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(VMH.comb, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(VMH.comb, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()