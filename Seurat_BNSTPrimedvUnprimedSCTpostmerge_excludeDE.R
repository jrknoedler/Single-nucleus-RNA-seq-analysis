#!/usr/bin/env Rscript
output <- "Seurat/BNST_POA_PrimedvUnprimedClusterTests/BNST_PrimedvUnprimed_excludeDE"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
library(matrixStats)
library(qlcMatrix)
library(igraph)
library(RANN)
PrimedBNST.data <- Read10X(data.dir= "MaleIntactVMH_retry/outs/filtered_feature_bc_matrix")
PrimedBNST <- CreateSeuratObject(counts = PrimedBNST.data, project = "PrimedBNST", min.cells=3, min.features=200)
UnprimedBNST.data <- Read10X(data.dir= "FemaleBNSTUnprimed/outs/filtered_feature_bc_matrix")
UnprimedBNST <- CreateSeuratObject(counts = UnprimedBNST.data, project = "UnprimedBNST", min.cells=3, min.features=200)
PrimedBNST$hormone <- "Primed"
UnprimedBNST$hormone <- "Unprimed"
mySeurat <- merge(PrimedBNST, y=c(UnprimedBNST), add.cell.ids=c("PrimedBNST","UnprimedBNST"), project="Merged")
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress="orig.ident")
genelist <- read.table("Genelists/BNST_PrimedvUnprimed.txt", header=FALSE)
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
mySeurat <- FindClusters(mySeurat, resolution=1.5)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_hormoneplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat,reduction="umap", label=TRUE, group.by="hormone")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))


