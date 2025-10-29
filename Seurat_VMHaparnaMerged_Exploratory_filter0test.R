#!/usr/bin/env Rscript

output <- "Seurat/VMH_PrimedvMale_AparnaNorm/VMH_MalevPrimed_aparnanorm_filter0test"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

Sex.integrated <- readRDS("Seurat/VMH_PrimedvMale_AparnaNorm/VMH_MalevPrimed_aparnanorm_filtered1.rds")
Sex.integrated <- subset(Sex.integrated, idents=c("0"), invert=TRUE)
genelist <- read.table("Genelists/VMH_MalevPrimedAll.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
head(Sex.integrated[[]])
hvg <- Sex.integrated@assays$RNA@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
Sex.integrated <- RunPCA(Sex.integrated, verbose=TRUE, features=hvg.final)
Sex.integrated <- RunUMAP(Sex.integrated, dims=1:30)
Sex.integrated <- FindNeighbors(Sex.integrated, dims=1:30)
Sex.integrated <- FindClusters(Sex.integrated, resolution=1.5)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(Sex.integrated, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexplit.pdf"))
DimPlot(Sex.integrated, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(Sex.integrated,reduction="umap", label=TRUE, group.by="sex")
dev.off()
Sex.integrated.markers <- FindAllMarkers(Sex.integrated, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(Sex.integrated.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- Sex.integrated.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Sex.integrated, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Sex.integrated, file=(paste0(output, ".rds")))
