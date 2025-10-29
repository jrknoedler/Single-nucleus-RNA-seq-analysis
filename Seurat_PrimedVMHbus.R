#!/usr/bin/env Rscript

sampleID <- "MaleVMH"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/PrimedVMH_Bus/PrimedVMH_Exploratory_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedVMH/combined_counts_em/", name="output", tcc=FALSE)
mySeurat <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=500)
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat
mySeurat <- SCTransform(mySeurat, verbose=TRUE)
mySeurat <- RunPCA(mySeurat)
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:10)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:10)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

