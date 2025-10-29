#!/usr/bin/env Rscript

sampleID <- "MaleVMH"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_Scale_DEPCA"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedVMH/combined_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=500)
Primed
Intact.data <- read_count_output("Kallisto_Bus_Output/IntactVMH/combined_counts_em/", name="output", tcc=FALSE)
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactMale", min.cells=3, min.features=500)
Intact
Unprimed.data <- read_count_output("Kallisto_Bus_Output/UnprimedVMH/combined_counts_em/", name="output", tcc=FALSE)
Unprimed <- CreateSeuratObject(counts = Unprimed.data, project = "UnprimedFemale", min.cells=3, min.features=500)
Unprimed
Intact$sex <- "Male"
Primed$sex <- "Female"
Unprimed$sex <- "Female"
Intact$status <- "IntactMale"
Primed$status <- "PrimedFemale"
Unprimed$status <- "UnprimedFemale"
Intact$hormone <- "Steroidal"
Primed$hormone <- "Steroidal"
Unprimed$hormone <- "Nonsteroidal"
mySeurat <- merge(Intact, y=c(Primed, Unprimed),add.cell.ids=c("IntactMale","PrimedFemale","UnprimedFemale"), project="Merged")
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat
genelist <- read.table("Genelists/VMH_allDEnew.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
DE <- read.table("Genelists/VMH_allDEnew.txt", header=FALSE)
DE <- unlist(DE)
DE <- as.matrix(DE)
mySeurat
head(mySeurat[[]])
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
mySeurat <- NormalizeData(mySeurat)
mySeurat <- FindVariableFeatures(mySeurat, selection.method="vst", nfeatures=2000)
mySeurat <- ScaleData(mySeurat)
mySeurat <- RunPCA(mySeurat, features=filtered.genelist)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_Cckarfeatureplot_umap.pdf"))
FeaturePlot(mySeurat, features="Cckar", reduction="umap", cols=c("yellow","purple"), label=FALSE)
dev.off()
pdf(paste0(output,"_UMAP_sexplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat,reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_statussplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="status")
dev.off()
pdf(paste0(output,"_UMAP_statuslabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="status")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

