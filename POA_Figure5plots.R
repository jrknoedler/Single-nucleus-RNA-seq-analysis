#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_CelltypesKiss1regressed_Fig5April"
regoutput <- "RegulonsCompiled/POA/Sexmarkersperclust/"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
pdf(file=paste0(output,"NewUMAPsize5.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=5)
dev.off()
pdf(file=paste0(output,"NewUMAPsize6.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=6)
dev.off()
pdf(file=paste0(output,"NewUMAPsize7.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=7)
dev.off()
pdf(file=paste0(output,"NewUMAPsize8.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=8)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)


pdf(file=paste0(output,"Markervln.pdf"), width=40, height=20)
VlnPlot(mySeurat, features=c("Esr1","Gad1","Slc17a6","Fbn2"), ncol=1, pt.size=0)
dev.off()
pdf(file=paste0(output,"Fbnftr.pdf"))
FeaturePlot(mySeurat, features=c("Fbn2"), cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output,"Fbnftr_hormonesplit.pdf"), width=30)
FeaturePlot(mySeurat, features=c("Fbn2"), split.by="Hormone", cols=c("light blue","red"))
dev.off()
pdf(file=paste0(output,"sDEGvln.pdf"), width=40, height=50)
VlnPlot(mySeurat, features=c("Col25a1","Pdzrn4","Nrip1","Shisa6","Ecel1","Esr2","Ryr3","Xkr4","Crim1","Tnxb","Prlr","Etl4"), ncol=1, pt.size=0)
dev.off()
male <- subset(mySeurat, idents=c("24","13","19","33"))
male.markers <- FindAllMarkers(male, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.30, test.use="MAST",  min.diff.pct = 0.2)
write.csv(male.markers, file=paste0(output,"nmalemarkers.csv"))


pdf(file=paste0(output,"brainstormvln_male.pdf"), width=10, height=100)
VlnPlot(male, features=c("Esr1","Ar","Gad1","Slc17a6","Slc18a2","Sst","Cartpt","Gal","Cck","Esr2","Tac1","Tac2","Th","Maob","Kiss1","Pappa","Col25a1","Npy","Fbn2","Calcr","Tacr1","Tacr3","Npy1r","Npy2r"), ncol=1, pt.size=0)
dev.off()


male <- ScaleData(male, features=rownames(male), do.center=TRUE, do.scale=TRUE)
pdf(paste0(output,"_Maletop10heatmap.pdf"))
top10 <- male.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(male, features=top10$gene, raster=FALSE) 
dev.off()

top10genes <- top10$gene
top10genes <- unique(top10genes)
Avg <- AverageExpression(male, assays=c("RNA"), features=top10genes)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/MakeMarkertop10_clusteravgs.csv")
pdf(paste0(output,"_Maletop5Markerheatmap.pdf"))
top5 <- male.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(male, features=top5$gene, raster=FALSE) 
dev.off()

top5genes <- top5$gene
top5genes <- unique(top5genes)
Avg <- AverageExpression(male, assays=c("RNA"), features=top5genes)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/MakeMarkertop5_clusteravgs.csv")