#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_Aromataseneurons"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_MalevUnprimed_Exploratory/BNST_MalevUnprimed_doubletsremoved_filtered2.rds")
mySeurat
filteredSeurat <- subset(mySeurat, idents=c("1","6","7","12","14","15","17","20"))
DefaultAssay(filteredSeurat) <- "integrated"
filteredSeurat <- RunUMAP(filteredSeurat, dims=1:25)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(filteredSeurat, reduction="umap", label=TRUE, group.by="sex")
dev.off()
DefaultAssay(filteredSeurat) <- "RNA"
NormalizeData(filteredSeurat)
pdf(paste0(output,"_Cyp19a1featureplot.pdf"))
FeaturePlot(mySeurat, blend.threshold=0.1, min.cutoff="q10",max.cutoff="q90", cols=c("lightgrey","red"), features=c("Cyp19a1"), split.by="sex")
dev.off()
clusters <-  c("1","6","7","12","14","15","17","20")
for (i in clusters){
try({
filteredSeurat.markers <- FindConservedMarkers(filteredSeurat, ident.1= i, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, grouping.var = "sex")
write.csv(filteredSeurat.markers, file=paste0(output,i,"_allposmarkers.csv"))
})
}
saveRDS(filteredSeurat, file=(paste0(output, ".rds")))
