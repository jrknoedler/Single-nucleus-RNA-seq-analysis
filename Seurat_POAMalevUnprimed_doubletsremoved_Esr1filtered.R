#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_MalevUnprimedMerged_Exploratory/POA_MalevUnprimed_Esr1filtered"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_MalevUnprimedMerged_Exploratory/MaleUnprimed_doubletsremoved_filtered5.rds")
mySeurat
filteredSeurat <- subset(mySeurat, idents=c("1","14","8","13","5","18","3","15","28","4","22"))
DefaultAssay(filteredSeurat) <- "integrated"
filteredSeurat <- RunUMAP(filteredSeurat, dims=1:30)
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
clusters <-  c("1","14","8","13","5","18","3","15","28","4","22")
for (i in clusters){
try({
filteredSeurat.markers <- FindConservedMarkers(filteredSeurat, ident.1= i, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, grouping.var = "sex")
write.csv(filteredSeurat.markers, file=paste0(output,i,"_allposmarkers.csv"))
})
}
saveRDS(filteredSeurat, file=(paste0(output, ".rds")))
