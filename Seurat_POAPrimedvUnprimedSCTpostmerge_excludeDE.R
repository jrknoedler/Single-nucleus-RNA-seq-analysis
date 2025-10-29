#!/usr/bin/env Rscript
output <- "Seurat/BNST_POA_PrimedvUnprimedClusterTests/POA_PrimedvUnprimed_ExcludedDE"

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
PrimedPOA.data <- Read10X(data.dir= "MaleMeAIntact/outs/filtered_feature_bc_matrix")
PrimedPOA <- CreateSeuratObject(counts = PrimedPOA.data, project = "PrimedPOA", min.cells=3, min.features=200)
UnprimedPOA.data <- Read10X(data.dir= "FemalePOAUnprimed/outs/filtered_feature_bc_matrix")
UnprimedPOA <- CreateSeuratObject(counts = UnprimedPOA.data, project = "UnprimedPOA", min.cells=3, min.features=200)
PrimedPOA$hormone <- "Primed"
UnprimedPOA$hormone <- "Unprimed"
mySeurat <- merge(PrimedPOA, y=c(UnprimedPOA), add.cell.ids=c("PrimedPOA","UnprimedPOA"), project="Merged")
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress="orig.ident")
genelist <- read.table("Genelists/POA_PrimedvUnprimed.txt", header=FALSE)
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
mySeurat <- RunPCA(mySeurat, features = hvg.final)
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.5)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
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


