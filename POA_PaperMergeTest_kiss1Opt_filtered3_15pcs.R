#!/usr/bin/env Rscript

output <- "Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL_15pcs_res0.8"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered2_30pcs_res1.2.rds")
mySeurat<- subset(mySeurat, idents=c("22","37"), invert=TRUE)
mySeurat <- subset(mySeurat, subset = nCount_RNA < 60000)
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
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:15)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:15)
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
pdf(paste0(output,"_Esr1Vln.pdf"), width=14)
VlnPlot(mySeurat, features=c("Esr1"), pt.size=0)
dev.off()
pdf(paste0(output,"UMI1Vln.pdf"), width=14)
VlnPlot(mySeurat, features=c("nCount_RNA"), pt.size=0)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(mySeurat, file=paste0(output, ".rds"))