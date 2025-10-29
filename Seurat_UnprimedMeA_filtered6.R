#!/usr/bin/env Rscript

output <- "Seurat/UnprimedMeA_Exploratory/UnprimedMeA_1000cutoff_filtered6_40pcs"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/UnprimedMeA_Exploratory/UnprimedMeA_1000cutoff_filtered5.rds")
mySeurat <- subset(mySeurat, idents=c("0"), invert=TRUE)
mySeurat
mySeurat <- RunPCA(mySeurat, feature=VariableFeatures(object=mySeurat))
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:40)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:40)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

