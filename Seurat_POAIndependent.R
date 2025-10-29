#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_IndependentRound1_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedPOA/combined_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=500)
Primed
Intact.data <- read_count_output("Kallisto_Bus_Output/IntactPOA/combined_counts_em/", name="output", tcc=FALSE)
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactMale", min.cells=3, min.features=500)
Intact
Unprimed.data <- read_count_output("Kallisto_Bus_Output/UnprimedPOA/combined_counts_em/", name="output", tcc=FALSE)
Unprimed <- CreateSeuratObject(counts = Unprimed.data, project = "UnprimedFemale", min.cells=3, min.features=500)
Unprimed
Castrate.data <- read_count_output("Kallisto_Bus_Output/CastratePOA/combined_counts_em/", name="output", tcc=FALSE)
Castrate <- CreateSeuratObject(counts = Castrate.data, project = "CastrateMale", min.cells=3, min.features=500)
Castrate

#cluster and save sample 1
Primed <- subset(Primed, subset=nCount_RNA < 60000)
Primed
Primed <- SCTransform(Primed, verbose=TRUE)
Primed <- RunPCA(Primed)
warnings()
Primed <- FindNeighbors(Primed, dims=1:30)
Primed <- FindClusters(Primed, resolution=1)
Primed <- RunUMAP(Primed, dims=1:30)
pdf(paste0(output,"_Primed_UMAP.pdf"))
DimPlot(Primed, reduction="umap", label=TRUE)
dev.off()
Primed.markers <- FindAllMarkers(Primed, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Primed.markers, file=paste0(output,"_Primed_allposmarkers.csv"))
pdf(paste0(output,"_Primed_TopMarkerheatmap.pdf"))
top10 <- Primed.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Primed, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Primed, file=(paste0(output, "_Primed.rds")))

#cluster and save sample 2
Unprimed <- subset(Unprimed, subset=nCount_RNA < 60000)
Unprimed <- SCTransform(Unprimed, verbose=TRUE)
Unprimed <- RunPCA(Unprimed)
warnings()
Unprimed <- FindNeighbors(Unprimed, dims=1:30)
Unprimed <- FindClusters(Unprimed, resolution=1)
Unprimed <- RunUMAP(Unprimed, dims=1:30)
pdf(paste0(output,"_Unprimed_UMAP.pdf"))
DimPlot(Unprimed, reduction="umap", label=TRUE)
dev.off()
Unprimed.markers <- FindAllMarkers(Unprimed, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Unprimed.markers, file=paste0(output,"_Unprimed_allposmarkers.csv"))
pdf(paste0(output,"_Unprimed_TopMarkerheatmap.pdf"))
top10 <- Unprimed.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Unprimed, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Unprimed, file=(paste0(output, "_Unprimed.rds")))

#cluster and save sample 3
Intact <- subset(Intact, subset=nCount_RNA < 60000)
Intact <- SCTransform(Intact, verbose=TRUE)
Intact <- RunPCA(Intact)
warnings()
Intact <- FindNeighbors(Intact, dims=1:30)
Intact <- FindClusters(Intact, resolution=1)
Intact <- RunUMAP(Intact, dims=1:30)
pdf(paste0(output,"_Intact_UMAP.pdf"))
DimPlot(Intact, reduction="umap", label=TRUE)
dev.off()
Intact.markers <- FindAllMarkers(Intact, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Intact.markers, file=paste0(output,"_Intact_allposmarkers.csv"))
pdf(paste0(output,"_Intact_TopMarkerheatmap.pdf"))
top10 <- Intact.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Intact, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Intact, file=(paste0(output, "_Intact.rds")))

#cluster and save sample 4
Castrate <- subset(Castrate, subset=nCount_RNA < 60000)
Castrate <- SCTransform(Castrate, verbose=TRUE)
Castrate <- RunPCA(Castrate)
warnings()
Castrate <- FindNeighbors(Castrate, dims=1:30)
Castrate <- FindClusters(Castrate, resolution=1)
Castrate <- RunUMAP(Castrate, dims=1:30)
pdf(paste0(output,"_Castrate_UMAP.pdf"))
DimPlot(Castrate, reduction="umap", label=TRUE)
dev.off()
Castrate.markers <- FindAllMarkers(Castrate, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Castrate.markers, file=paste0(output,"_Castrate_allposmarkers.csv"))
pdf(paste0(output,"_Castrate_TopMarkerheatmap.pdf"))
top10 <- Castrate.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Castrate, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Castrate, file=(paste0(output, "_Castrate.rds")))

