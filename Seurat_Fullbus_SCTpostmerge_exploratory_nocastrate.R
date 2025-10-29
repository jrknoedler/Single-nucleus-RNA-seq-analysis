#!/usr/bin/env Rscript

sampleID <- "MaleVMH"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/FULLBUS/FullBusMerge_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
PrimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1__Primed.rds")
head(PrimedVMH[[]])
PrimedVMH$orig.ident <- NULL
head(PrimedVMH[[]])
IntactVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1__Intact.rds")
IntactVMH$orig.ident <- NULL
UnprimedVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1__Unprimed.rds")
UnprimedVMH$orig.ident <- NULL
PrimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1__Primed.rds")
PrimedBNST$orig.ident <- NULL
IntactBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1__Intact.rds")
IntactBNST$orig.ident <- NULL
UnprimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1__Unprimed.rds")
UnprimedBNST$orig.ident <- NULL
PrimedPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentFiltered1__Primed.rds")
PrimedPOA$orig.ident <- NULL
IntactPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentFiltered1__Intact.rds")
IntactPOA$orig.ident <- NULL
UnprimedPOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentFiltered1__Unprimed.rds")
UnprimedPOA$orig.ident <- NULL
CastratePOA <- readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentFiltered1__Castrate.rds")
CastratePOA$orig.ident <- NULL
PrimedMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered1__Primed.rds")
PrimedMeA$orig.ident <- NULL
IntactMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered1__Intact.rds")
IntactMeA$orig.ident <- NULL
UnprimedMeA <- readRDS("Seurat/UnprimedMeA_Exploratory/UnprimedMeA_1000cutoff_filtered3.rds")
UnprimedMeA$orig.ident <- NULL
CastrateMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered1__Castrate.rds")
CastrateMeA$orig.ident <- NULL
CastrateVMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1__Castrate_filtered2.rds")
CastrateVMH$orig.ident <- NULL
IntactVMH$region <- "VMH"
IntactBNST$region <- "BNST"
IntactPOA$region <- "POA"
IntactMeA$region <- "MeA"
PrimedVMH$region <- "VMH"
PrimedBNST$region <- "BNST"
PrimedPOA$region <- "POA"
PrimedMeA$region <- "MeA"
UnprimedVMH$region <- "VMH"
UnprimedBNST$region <- "BNST"
UnprimedPOA$region <- "POA"
UnprimedMeA$region <- "MeA"
CastrateVMH$region <- "VMH"
CastratePOA$region <- "POA"
CastrateMeA$region <- "MeA"

IntactVMH$status <- "IntactVMH"
IntactBNST$status <- "IntactBNST"
IntactPOA$status <- "IntactPOA"
IntactMeA$status <- "IntactMeA"
PrimedVMH$status <- "PrimedVMH"
PrimedBNST$status <- "PrimedBNST"
PrimedPOA$status <- "PrimedPOA"
PrimedMeA$status <- "PrimedMeA"
UnprimedVMH$status <- "UnprimedVMH"
UnprimedBNST$status <- "UnprimedBNST"
UnprimedPOA$status <- "UnprimedPOA"
UnprimedMeA$status <- "UnprimedMeA"
CastrateVMH$status <- "CastrateVMH"
CastratePOA$status <- "CastratePOA"
CastrateMeA$status <- "CastrateMeA"

Shared.features <- SelectIntegrationFeatures(object.list=list(IntactVMH,PrimedVMH, UnprimedVMH, CastrateVMH, IntactBNST,PrimedBNST,UnprimedBNST,IntactPOA,PrimedPOA,UnprimedPOA,CastratePOA,IntactMeA,PrimedMeA,UnprimedMeA,CastrateMeA), nfeatures=3000)
Shared.list <- PrepSCTIntegration(object.list=list(IntactVMH,PrimedVMH, UnprimedVMH, CastrateVMH, IntactBNST,PrimedBNST,UnprimedBNST,IntactPOA,PrimedPOA,UnprimedPOA,CastratePOA,IntactMeA,PrimedMeA,UnprimedMeA,CastrateMeA), anchor.features=Shared.features, verbose=TRUE)
Shared.anchors <- FindIntegrationAnchors(object.list=Shared.list, normalization.method="SCT", anchor.features=Shared.features, verbose=TRUE)
mySeurat <- IntegrateData(anchorset=Shared.anchors, normalization.method="SCT", verbose=FALSE)
mySeurat <- RunPCA(mySeurat, npcs=100)
mySeurat
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:100)
mySeurat <- FindClusters(mySeurat, resolution=2)
mySeurat <- RunUMAP(mySeurat, dims=1:100)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_regionplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="region")
dev.off()
pdf(paste0(output,"_UMAP_regionlabel.pdf"))
DimPlot(mySeurat,reduction="umap", label=TRUE, group.by="region")
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output,"_Slc17a6feature.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Slc17a6", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Esr1feature.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Esr1", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Slc32a1feature.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Slc32a1", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Olig1feature.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Olig1", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Mobpfeature.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Mobp", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Gfapfeature.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Gfap", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Slc1a3.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Slc1a3", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Cx3cr1.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Cx3cr1", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Inpp5d.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Inpp5d", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Syn1.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Syn1", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Syn2.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Syn2", cols=c("light blue", "red"))
dev.off()
pdf(paste0(output,"_Syn3.pdf")
FeaturePlot(mySeurat, reduction="umap", features="Syn3", cols=c("light blue", "red"))
dev.off()
mySeurat <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

