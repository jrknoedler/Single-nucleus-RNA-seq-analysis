#!/usr/bin/env Rscript

output <- "Seurat/MeA_IndependentAnalysis/MeA_MergedCCA_Filtered2_res0.8_20pcs"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
#library(DoubletFinder)

mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_MergedCCA_Filtered1_res1.2.rds")
mySeurat <- subset(mySeurat, idents=c("30"), invert=TRUE)
#mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
#mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.Malat1"))
DefaultAssay(mySeurat) <- "integrated"
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$integrated@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
mySeurat <- RunPCA(mySeurat, features=hvg.final)
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:20)
mySeurat <- FindClusters(mySeurat, resolution=0.8)
mySeurat <- RunUMAP(mySeurat, dims=1:20)
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
pdf(paste0(output,"_UMIvln.pdf"), width=12)
VlnPlot(mySeurat,features=c("nCount_RNA"), pt.size=1)
dev.off() 
pdf(paste0(output,"_Esr1vln.pdf"), width=12)
VlnPlot(mySeurat,features=c("Esr1"), pt.size=0)
dev.off() 
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

