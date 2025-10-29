#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/WinterBreakMergetest/Aparnamerge_100pcs"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
library(matrixStats)
library(qlcMatrix)
library(igraph)
library(RANN)
MaleBNST.data <- Read10X(data.dir = "MaleBNST_Runscombined/outs/filtered_feature_bc_matrix")
MaleBNST <- CreateSeuratObject(counts = MaleBNST.data, project = "MaleBNST", min.cells=3, min.features=200)
MaleBNSTtotal <- colSums(MaleBNST[["RNA"]]@counts[rownames(MaleBNST[["RNA"]]@counts),])
MaleBNST_median <- median(MaleBNSTtotal)
MaleBNSTmaxes <- colMax(MaleBNST[["RNA"]]@counts[rownames(MaleBNST[["RNA"]]@counts),])
MaleBNSTnorm <- MaleBNST[["RNA"]]@counts[rownames(MaleBNST[["RNA"]]@counts),]/as.numeric(MaleBNSTmaxes)
MaleBNSTcorrected <- MaleBNSTnorm * MaleBNST_median
MaleBNSTfinal <- CreateSeuratObject(counts = MaleBNSTcorrected, project = "IntactMaleBNST")
head(MaleBNSTfinal[[]])
MalePOA.data <- Read10X(data.dir= "MalePOAIntact/outs/filtered_feature_bc_matrix")
MalePOA <- CreateSeuratObject(counts = MalePOA.data, project = "MalePOA", min.cells=3, min.features=200)
MalePOAtotal <- colSums(MalePOA[["RNA"]]@counts[rownames(MalePOA[["RNA"]]@counts),])
MalePOA_median <- median(MalePOAtotal)
MalePOAmaxes <- colMax(MalePOA[["RNA"]]@counts[rownames(MalePOA[["RNA"]]@counts),])
MalePOAnorm <- MalePOA[["RNA"]]@counts[rownames(MalePOA[["RNA"]]@counts),]/as.numeric(MalePOAmaxes)
MalePOAcorrected <- MalePOAnorm * MalePOA_median
MalePOAfinal <- CreateSeuratObject(counts = MalePOAcorrected, project = "IntactMalePOA")
MalePOAfinal
PrimedBNST.data <- Read10X(data.dir= "MaleIntactVMH_retry/outs/filtered_feature_bc_matrix")
PrimedBNST <- CreateSeuratObject(counts = PrimedBNST.data, project = "PrimedBNST", min.cells=3, min.features=200)
PrimedBNSTtotal <- colSums(PrimedBNST[["RNA"]]@counts[rownames(PrimedBNST[["RNA"]]@counts),])
PrimedBNST_median <- median(PrimedBNSTtotal)
PrimedBNSTmaxes <- colMax(PrimedBNST[["RNA"]]@counts[rownames(PrimedBNST[["RNA"]]@counts),])
PrimedBNSTnorm <- PrimedBNST[["RNA"]]@counts[rownames(PrimedBNST[["RNA"]]@counts),]/as.numeric(PrimedBNSTmaxes)
PrimedBNSTcorrected <- PrimedBNSTnorm * PrimedBNST_median
PrimedBNSTfinal <- CreateSeuratObject(counts = PrimedBNSTcorrected, project = "PrimedBNST")
PrimedBNSTfinal
PrimedPOA.data <- Read10X(data.dir= "MaleMeAIntact/outs/filtered_feature_bc_matrix")
PrimedPOA <- CreateSeuratObject(counts = PrimedPOA.data, project = "PrimedPOA", min.cells=3, min.features=200)
PrimedPOAtotal <- colSums(PrimedPOA[["RNA"]]@counts[rownames(PrimedPOA[["RNA"]]@counts),])
PrimedPOA_median <- median(PrimedPOAtotal)
PrimedPOAmaxes <- colMax(PrimedPOA[["RNA"]]@counts[rownames(PrimedPOA[["RNA"]]@counts),])
PrimedPOAnorm <- PrimedPOA[["RNA"]]@counts[rownames(PrimedPOA[["RNA"]]@counts),]/as.numeric(PrimedPOAmaxes)
PrimedPOAcorrected <- PrimedPOAnorm * PrimedPOA_median
PrimedPOAfinal <- CreateSeuratObject(counts = PrimedPOAcorrected, project = "PrimedPOA")
PrimedPOAfinal
MaleMeA.data <- Read10X(data.dir= "FemaleBNSTPrimed_Putative_retry/outs/filtered_feature_bc_matrix")
MaleMeA <- CreateSeuratObject(counts = MaleMeA.data, project = "MaleMeA", min.cells=3, min.features=200)
MaleMeAtotal <- colSums(MaleMeA[["RNA"]]@counts[rownames(MaleMeA[["RNA"]]@counts),])
MaleMeA_median <- median(MaleMeAtotal)
MaleMeAmaxes <- colMax(MaleMeA[["RNA"]]@counts[rownames(MaleMeA[["RNA"]]@counts),])
MaleMeAnorm <- MaleMeA[["RNA"]]@counts[rownames(MaleMeA[["RNA"]]@counts),]/as.numeric(MaleMeAmaxes)
MaleMeAcorrected <- MaleMeAnorm * MaleMeA_median
MaleMeAfinal <- CreateSeuratObject(counts = MaleMeAcorrected, project = "MaleMeA")
MaleMeAfinal
MaleVMH.data <- Read10X(data.dir = "FemalePOAPrimed_Putative2/outs/filtered_feature_bc_matrix")
MaleVMH <- CreateSeuratObject(counts = MaleVMH.data, project = "MaleVMH", min.cells=3, min.features=200)
MaleVMHtotal <- colSums(MaleVMH[["RNA"]]@counts[rownames(MaleVMH[["RNA"]]@counts),])
MaleVMH_median <- median(MaleVMHtotal)
MaleVMHmaxes <- colMax(MaleVMH[["RNA"]]@counts[rownames(MaleVMH[["RNA"]]@counts),])
MaleVMHnorm <- MaleVMH[["RNA"]]@counts[rownames(MaleVMH[["RNA"]]@counts),]/as.numeric(MaleVMHmaxes)
MaleVMHcorrected <- MaleVMHnorm * MaleVMH_median
MaleVMHfinal <- CreateSeuratObject(counts = MaleVMHcorrected, project = "MaleVMH")
MaleVMHfinal
PrimedMeA.data <- Read10X(data.dir= "FemaleMeAPrimed/outs/filtered_feature_bc_matrix")
PrimedMeA <- CreateSeuratObject(counts = PrimedMeA.data, project = "PrimedMeA", min.cells=3, min.features=200)
PrimedMeAtotal <- colSums(PrimedMeA[["RNA"]]@counts[rownames(PrimedMeA[["RNA"]]@counts),])
PrimedMeA_median <- median(PrimedMeAtotal)
PrimedMeAmaxes <- colMax(PrimedMeA[["RNA"]]@counts[rownames(PrimedMeA[["RNA"]]@counts),])
PrimedMeAnorm <- PrimedMeA[["RNA"]]@counts[rownames(PrimedMeA[["RNA"]]@counts),]/as.numeric(PrimedMeAmaxes)
PrimedMeAcorrected <- PrimedMeAnorm * PrimedMeA_median
PrimedMeAfinal <- CreateSeuratObject(counts = PrimedMeAcorrected, project = "PrimedMeA")
PrimedMeAfinal
PrimedVMH.data <- Read10X(data.dir= "FemaleVMHPrimed/outs/filtered_feature_bc_matrix")
PrimedVMH <- CreateSeuratObject(counts = PrimedVMH.data, project = "PrimedVMH", min.cells=3, min.features=200)
PrimedVMHtotal <- colSums(PrimedVMH[["RNA"]]@counts[rownames(PrimedVMH[["RNA"]]@counts),])
PrimedVMH_median <- median(PrimedVMHtotal)
PrimedVMHmaxes <- colMax(PrimedVMH[["RNA"]]@counts[rownames(PrimedVMH[["RNA"]]@counts),])
PrimedVMHnorm <- PrimedVMH[["RNA"]]@counts[rownames(PrimedVMH[["RNA"]]@counts),]/as.numeric(PrimedVMHmaxes)
PrimedVMHcorrected <- PrimedVMHnorm * PrimedVMH_median
PrimedVMHfinal <- CreateSeuratObject(counts = PrimedVMHcorrected, project = "PrimedVMH")
PrimedVMHfinal
UnprimedBNST.data <- Read10X(data.dir= "FemaleBNSTUnprimed/outs/filtered_feature_bc_matrix")
UnprimedBNST <- CreateSeuratObject(counts = UnprimedBNST.data, project = "UnprimedBNST", min.cells=3, min.features=200)
UnprimedBNSTtotal <- colSums(UnprimedBNST[["RNA"]]@counts[rownames(UnprimedBNST[["RNA"]]@counts),])
UnprimedBNST_median <- median(UnprimedBNSTtotal)
UnprimedBNSTmaxes <- colMax(UnprimedBNST[["RNA"]]@counts[rownames(UnprimedBNST[["RNA"]]@counts),])
UnprimedBNSTnorm <- UnprimedBNST[["RNA"]]@counts[rownames(UnprimedBNST[["RNA"]]@counts),]/as.numeric(UnprimedBNSTmaxes)
UnprimedBNSTcorrected <- UnprimedBNSTnorm * UnprimedBNST_median
UnprimedBNSTfinal <- CreateSeuratObject(counts = UnprimedBNSTcorrected, project = "UnprimedBNST")
UnprimedBNSTfinal
UnprimedPOA.data <- Read10X(data.dir= "FemalePOAUnprimed/outs/filtered_feature_bc_matrix")
UnprimedPOA <- CreateSeuratObject(counts = UnprimedPOA.data, project = "UnprimedPOA", min.cells=3, min.features=200)
UnprimedPOAtotal <- colSums(UnprimedPOA[["RNA"]]@counts[rownames(UnprimedPOA[["RNA"]]@counts),])
UnprimedPOA_median <- median(UnprimedPOAtotal)
UnprimedPOAmaxes <- colMax(UnprimedPOA[["RNA"]]@counts[rownames(UnprimedPOA[["RNA"]]@counts),])
UnprimedPOAnorm <- UnprimedPOA[["RNA"]]@counts[rownames(UnprimedPOA[["RNA"]]@counts),]/as.numeric(UnprimedPOAmaxes)
UnprimedPOAcorrected <- UnprimedPOAnorm * UnprimedPOA_median
UnprimedPOAfinal <- CreateSeuratObject(counts = UnprimedPOAcorrected, project = "UnprimedPOA")
UnprimedPOAfinal
MaleVMHfinal$Sex <- "Male"
MaleVMHfinal$Region <- "VMH"
MaleVMHfinal$Batch <- "3"
MaleVMHfinal$Hormone <- "Intact"
MaleBNSTfinal$Sex <- "Male"
MaleBNSTfinal$Region <- "BNST"
MaleBNSTfinal$Batch <- "1"
MaleBNSTfinal$Hormone <- "Intact"
MalePOAfinal$Sex <- "Male"
MalePOAfinal$Region <- "POA"
MalePOAfinal$Batch <- "2"
MalePOAfinal$Hormone <- "Intact"
MaleMeAfinal$Sex <- "Male"
MaleMeAfinal$Region <- "MeA"
MaleMeAfinal$Batch <- "3"
MaleMeAfinal$Hormone <- "Intact"
PrimedBNSTfinal$Sex <- "Female"
PrimedBNSTfinal$Region <- "BNST"
PrimedBNSTfinal$Batch <- "3"
PrimedBNSTfinal$Hormone <- "Primed"
PrimedPOAfinal$Sex <- "Female"
PrimedPOAfinal$Region <- "POA"
PrimedPOAfinal$Batch <- "3"
PrimedPOAfinal$Hormone <- "Primed"
PrimedMeAfinal$Sex <- "Female"
PrimedMeAfinal$Region <- "MeA"
PrimedMeAfinal$Batch <- "3"
PrimedMeAfinal$Hormone <- "Primed"
PrimedVMHfinal$Sex <- "Female"
PrimedVMHfinal$Region <- "VMH"
PrimedVMHfinal$Batch <- "3"
PrimedVMHfinal$Hormone <- "Primed"
UnprimedPOAfinal$Sex <- "Female"
UnprimedPOAfinal$Region <- "POA"
UnprimedPOAfinal$Batch <- "2"
UnprimedPOAfinal$Hormone <- "Unprimed"
UnprimedBNSTfinal$Sex <- "Female"
UnprimedBNSTfinal$Region <- "BNST"
UnprimedBNSTfinal$Batch <- "2"
UnprimedBNSTfinal$Hormone <- "Unprimed"
mySeurat <- merge(MaleBNSTfinal, y=c(MalePOAfinal, MaleMeAfinal, MaleVMHfinal, UnprimedPOAfinal, UnprimedBNSTfinal, PrimedPOAfinal, PrimedBNSTfinal, PrimedMeAfinal, PrimedVMHfinal), add.cell.ids=c("MaleBNSTfinal","MalePOAfinal", "MaleMeAfinal","MaleVMHfinal","PrimedBNSTfinal","PrimedPOAfinal","PrimedMeAfinal","PrimedVMHfinal","UnprimedBNSTfinal","UnprimedPOAfinal"), project="FullHouseMerge")
mySeurat
head(mySeurat[[]])
mySeurat <- NormalizeData(mySeurat)
mySeurat <- FindVariableFeatures(mySeurat, selection.method="vst", nfeatures=2000)
mySeurat <- ScaleData(mySeurat, vars.to.regress="orig.ident")
mySeurat <- RunPCA(mySeurat, npcs=100)
pdf(paste0(output,"_PCAplot.pdf"))
DimPlot(mySeurat, reduction="pca", group.by="orig.ident")
dev.off()
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:100)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:100)
mySeurat <- FindClusters(mySeurat, resolution=1.5)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_origsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="orig.ident")
dev.off()
pdf(paste0(output,"_UMAP_origlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="orig.ident")
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_regionsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Region")
dev.off()
pdf(paste0(output,"_UMAP_regionlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Region")
dev.off()
pdf(paste0(output,"_UMAP_batchsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Batch")
dev.off()
pdf(paste0(output,"_UMAP_batchlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Batch")
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
FeaturePlot(mySeurat, features=c("Xist"), split.by="orig.ident")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
