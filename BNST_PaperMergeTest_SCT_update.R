#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_SCT20pcsres0.8sanitycheck"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
PrimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Prepaperrecluster_2_Primed.rds")
PrimedBNST <- subset(PrimedBNST, idents=c("2"), invert=TRUE)
PrimedBNST
MaleBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Prepaperrecluster_2_Intact.rds")
MaleBNST <- subset(MaleBNST, idents=c("25","0"), invert=TRUE)
MaleBNST
UnprimedBNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Prepaperrecluster_2_Unprimed.rds")
UnprimedBNST <- subset(UnprimedBNST, idents=c("0"), invert=TRUE)
UnprimedBNST
#PrimedBNST <- SCTransform(PrimedBNST)
#MaleBNST <- SCTransform(MaleBNST)
#UnprimedBNST <- SCTransform(UnprimedBNST)
MaleBNST$sex <- "Male"
MaleBNST$Hormone <- "Intact"
PrimedBNST$sex <- "Female"
PrimedBNST$Hormone <- "Primed"
UnprimedBNST$sex <- "Female"
UnprimedBNST$Hormone <- "Unprimed"
list <- c(MaleBNST,UnprimedBNST,PrimedBNST)
BNST.features <- SelectIntegrationFeatures(object.list=list, nfeatures=3000)
BNST.list <- PrepSCTIntegration(list, anchor.features=BNST.features, verbose=TRUE)
BNST.anchors <- FindIntegrationAnchors(BNST.list, normalization.method="SCT", anchor.features=BNST.features, verbose=TRUE)
mySeurat <- IntegrateData(anchorset=BNST.anchors, normalization.method="SCT", verbose=FALSE)
mySeurat <- RunPCA(mySeurat)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:20)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:20)
mySeurat <- FindClusters(mySeurat, resolution=0.8)
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
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(file=paste0(output,"Marker_Vln.pdf"), width=40, height=15)
VlnPlot(mySeurat, features=c("Gad1","Slc17a7","Esr1","Cyp19a1","Tac1"), pt.size=0)
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.20)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))