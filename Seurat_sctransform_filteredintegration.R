#!/usr/bin/env Rscript

output <- "Seurat/POA_MalevPrimed_Filteredintegration/POA_MalevPrimed_filtered"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
MaleBNST.data <- Read10X(data.dir = "MaleBNST_Runscombined/outs/filtered_feature_bc_matrix")
UnprimedBNST.data <- Read10X(data.dir = "FemaleBNSTUnprimed/outs/filtered_feature_bc_matrix")
MalePOA.data <- Read10X(data.dir = "MalePOAIntact/outs/filtered_feature_bc_matrix")
PrimedPOA.data <- Read10X(data.dir = "FemalePOAPrimed/outs/filtered_feature_bc_matrix")
UnprimedPOA.data <- Read10X(data.dir = "FemalePOAPrimed/outs/filtered_feature_bc_matrix")
MaleBNST <- CreateSeuratObject(counts = MaleBNST.data, project ="MaleBNST", min.cells=3, min.features=200)
UnprimedBNST <- CreateSeuratObject(counts = UnprimedBNST.data, project="UnprimedBNST", min.cells=3, min.features=200)
MalePOA <- CreateSeuratObject(counts = MalePOA.data, project="MalePOA", min.cells=3, min.features=200)
PrimedPOA <- CreateSeuratObject(counts = PrimedPOA.data, project="UnprimedPOA", min.cells=3, min.features=200)
UnprimedPOA <- CreateSeuratObject(counts = UnprimedPOA.data, project="UnprimedPOA", min.cells=3, min.features=200)
MaleBNST <- SCTransform(MaleBNST, verbose=TRUE)
UnprimedBNST <- SCTransform(UnprimedBNST, verbose=TRUE)
MalePOA <- SCTransform(MalePOA, verbose=TRUE)
PrimedPOA <- SCTransform(PrimedPOA, verbose=TRUE)
UnprimedPOA <- SCTransform(UnprimedPOA, verbose=TRUE)
MaleBNST$library <- "MaleBNST"
UnprimedBNST$library <- "UnprimedBNST"
MalePOA$library <- "MalePOA"
PrimedPOA$library <- "PrimedPOA"
UnprimedPOA$library <- "UnprimedPOA"
MaleBNST$sex <- "Male"
UnprimedBNST$sex <- "Female"
MalePOA$sex <- "Male"
PrimedPOA$sex <- "Female"
UnprimedPOA$sex <- "Female"
MaleBNST$region <- "BNST"
UnprimedBNST$region <- "BNST"
MalePOA$region <- "POA"
PrimedPOA$region <- "POA"
UnprimedPOA$region <- "POA"
Shared.features <- SelectIntegrationFeatures(object.list=list(MaleBNST,UnprimedBNST,MalePOA,PrimedPOA,UnprimedPOA), nfeatures=3000)
Shared.list <- PrepSCTIntegration(object.list=list(MaleBNST,UnprimedBNST,MalePOA,PrimedPOA,UnprimedPOA), anchor.features=Shared.features, verbose=TRUE)
Shared.anchors <- FindIntegrationAnchors(object.list=Shared.list, normalization.method="SCT", anchor.features=Shared.features, verbose=TRUE)
Full.integrated <- IntegrateData(anchorset=Shared.anchors, normalization.method="SCT", verbose=FALSE)
Full.integrated <- RunPCA(Full.integrated, verbose=TRUE)
Full.integrated <- RunUMAP(Full.integrated, dims=1:30)
Full.integrated <- FindNeighbors(Full.integrated, dims=1:30)
Full.integrated <- FindClusters(Full.integrated, resolution=1)
pdf("Seurat/Seurat_10x1stpast_FullUmap_Sctransform.pdf")
DimPlot(Full.integrated, reduction="umap", label=TRUE)
dev.off()
pdf("Seurat/Seurat_10x1stpast_FullUmap_Sctransform_librarysplit.pdf")
DimPlot(Full.integrated, reduction="umap", split.by="library")
dev.off()
pdf("Seurat/Seurat_10x1stpast_FullUmap_Sctransform_sexsplit.pdf")
DimPlot(Full.integrated, reduction="umap", split.by="sex")
dev.off()
pdf("Seurat/Seurat_10x1stpast_FullUmap_Sctransform_regionsplit.pdf")
DimPlot(Full.integrated, reduction="umap", split.by="region")
dev.off()
pdf("Seurat/Seurat_10x1stpast_FullUmap_Sctransform_regionslabel.pdf")
DimPlot(Full.integrated, reduction="umap", group.by="region")
dev.off()
pdf("Seurat/Seurat_10x1stpast_FullUmap_Sctransform_sexlabel.pdf")
DimPlot(Full.integrated, reduction="umap", group.by="sex")
dev.off()
pdf("Seurat/Seurat_10x1stpast_FullUmap_Sctransform_liblabel.pdf")
DimPlot(Full.integrated, reduction="umap", group.by="library")
dev.off()
cluster.markers <- FindAllMarkers(Full.integrated)
write.csv(cluster.markers, file="Seurat/MaleEsr1_Sctransformexploratory/MaleBNST_FullLiger_clusters_integrated.csv")
