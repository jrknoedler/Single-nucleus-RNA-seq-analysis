#!/usr/bin/env Rscript
output <- "Seurat/PrimedVMH_Bus/PrimedVMH_Exploratory_filtered4"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/PrimedVMH_Bus/PrimedVMH_Exploratory_filtered3.rds")
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, cols=c("Paired"))
dev.off()
pdf(paste0(output,"_TSNE.pdf"))
DimPlot(mySeurat, reduction="tsne", label=TRUE, cols=c("Paired"))
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(file=paste0(output, "Vln3.pdf"), width=14, height=28)
VlnPlot(mySeurat, features=c("Esr1","Onecut2","Klhl14","Sst","Nr5a1","Fezf1","Lncenc1","Abtb2","Coro6","Gal","Cckar"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output,"_Cckarfeatureplot_umap.pdf"))
FeaturePlot(mySeurat, features="Cckar", reduction="umap", cols=c("light blue","red"), label=FALSE)
dev.off()
pdf(paste0(output,"_Cckarfeatureplot_tsne.pdf"))
FeaturePlot(mySeurat, features="Cckar", reduction="tsne", cols=c("light blue","red"),  label=FALSE)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))

saveRDS(mySeurat, file=(paste0(output, ".rds")))

