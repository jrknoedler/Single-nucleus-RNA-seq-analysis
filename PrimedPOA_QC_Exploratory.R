#!/usr/bin/env Rscript

sampleID <- "PregPOA"
output <- "Seurat/Dec2019_MarkerTests/PrimedPOA_Test"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
Seurat.data <- read_count_output("Updated_Ref_Kb_Out/POA_Primed/output_mtx_em", name="output", tcc=FALSE)
mySeurat <- CreateSeuratObject(counts = Seurat.data, project = sampleID, min.cells=3, min.features=500)
mySeurat
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat <- SCTransform(mySeurat, verbose=TRUE)
mySeurat <- RunPCA(mySeurat)
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.5)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=60)
VlnPlot(mySeurat, features=c("Gad1","Slc17a6","Slc17a7","Nr5a1","Opn5","Cyp19a1","Cckar"), ncol=1, pt.size=0)
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))
