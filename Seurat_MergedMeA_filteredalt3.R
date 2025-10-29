#!/usr/bin/env Rscript

output <- "Seurat/MeA_IndependentAnalysis/MeA_MergedFinal_Filtered3alt"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
#library(DoubletFinder)

mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_MergedFinal_Filtered2alt.rds")
mySeurat <- subset(mySeurat, idents=c("4"), invert=TRUE)
mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.Malat1"))

genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final

mySeurat <- RunPCA(mySeurat, features=hvg.final)
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
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
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMIvln.pdf"), width=12)
VlnPlot(mySeurat,features=c("nCount_RNA"), pt.size=1)
dev.off() 
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#sweep.res.list <- paramSweep_v3(mySeurat, PCs=1:30, sct=TRUE)
#sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
#bcmvn <- find.pK(sweep.stats)
#bcmvn
#nExp_poi <- round(0.15*nrow(mySeurat@meta.data))
#mySeurat <- doubletFinder_v3(mySeurat, PCs=1:30, pN=0.25, pK=0.14, nExp=nExp_poi, reuse.pANN=FALSE, sct=TRUE)
#head(mySeurat[[]])
saveRDS(mySeurat, file=(paste0(output, ".rds")))
#pdf(file=paste0(output,"doublettest.pdf"))
#DimPlot(mySeurat, group.by="DF.classifications_0.25_0.14_1476")
#dev.off()


