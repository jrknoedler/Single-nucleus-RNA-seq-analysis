#!/usr/bin/env Rscript



output <- "Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_0.1CPM_hicat"
library(liger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(cowplot)
library(Matrix)
library(patchwork)
library(dplyr)
library(tidyverse)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")
Scratch.clusters <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Scratchcl.csv", header=TRUE)
scratch <- Scratch.clusters[["x"]]
unlist(scratch)
mySeurat <- AddMetaData(mySeurat, metadata=scratch, col.name="Scrattch_Clusters")

Idents(mySeurat) <- 'Scrattch_Clusters'
mySeurat@active.ident
pdf(file="Seurat/BNST_IndependentAnalysis/LabeledUMAPscratch.pdf")
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
Markers <- FindAllMarkers(mySeurat, test.use="MAST", min.pct=0.2, only.pos=TRUE)
write.csv(Markers, file="Seurat/BNST_IndependentAnalysis/BNST_Scratchclustermarkers.csv")

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/BNST_MvP_latest.txt")
unlist(Dimorphic)
Dimorphic <- as.matrix(Dimorphic)
Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)
mySeurat <- NormalizeData(mySeurat, scale.factor=1e6)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)
avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=Dimorphic.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in Dimorphic.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "Clusters_Per_SDEG.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "SDEGS_Per_Cluster.csv"))