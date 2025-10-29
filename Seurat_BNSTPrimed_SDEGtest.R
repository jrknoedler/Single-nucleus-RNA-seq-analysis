#!/usr/bin/env Rscript
output <- "Seurat/BNST_IndependentAnalysis/PrimedSDEG_Plottest"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1__Primed.rds")
head(mySeurat[[]])
genelist <- read.table("Genelists/BNST_MvP_Updated.txt")
unlist(genelist)
genelist <- as.matrix(genelist)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
SFARIgenes <- read.table("Genelists/SFARI.txt")
unlist(SFARIgenes)
SFARIgenes <- as.matrix(SFARIgenes)
filtered.SFARI <- intersect(SFARIgenes, genes.10x) 
Upmale <- read.table("Genelists/BNST_MvP_Upmale_latest.txt")
unlist(Upmale)
Upmale <- as.matrix(Upmale)
filtered.Upmale <- intersect(Upmale, genes.10x)
Upfemale <- read.table("Genelists/BNST_MvP_Upfemale_latest.txt")
unlist(Upfemale)
Upfemale <- as.matrix(Upfemale)
filtered.Upfemale <- intersect(Upfemale, genes.10x)
mySeurat[["percent.SDEG"]] <- PercentageFeatureSet(object=mySeurat, features=filtered.genelist)
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output,"_UMAP_pct_SDEG.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="percent.SDEG", cols = c("light blue", "red"))
dev.off()
pdf(paste0(output,"_UMAP_Esr1.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Esr1", cols = c("light blue", "red"))
dev.off()
pdf(paste0(output,"_UMAP_Cyp19a1.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Cyp19a1", cols = c("light blue", "red"))
dev.off()
pdf(paste0(output,"_UMAP_Ar.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Ar", cols = c("light blue", "red"))
dev.off()
mySeurat <- AddModuleScore(
object = mySeurat,
features = list(filtered.genelist),
name = 'SDEG.score'
)
mySeurat <- AddModuleScore(
object = mySeurat,
features = list(filtered.SFARI),
name = 'SFARI.score'
)
mySeurat <- AddModuleScore(
object = mySeurat,
features = list(filtered.Upmale),
name = 'Upmale.score'
)
mySeurat <- AddModuleScore(
object = mySeurat,
features = list(filtered.Upfemale),
name = 'Upfemale.score'
)
head(mySeurat[[]])
pdf(paste0(output,"_UMAP_SDEGmodule.pdf"))
FeaturePlot(mySeurat, features="SDEG.score1", cols = c("light blue", "red"))
dev.off()
pdf(paste0(output,"_UMAP_SFARImodule.pdf"))
FeaturePlot(mySeurat, features="SFARI.score1", cols = c("light blue", "red"))
dev.off()
pdf(paste0(output,"_UMAP_Upmalemodule.pdf"))
FeaturePlot(mySeurat, features="Upmale.score1", cols = c("light blue", "red"))
dev.off()
pdf(paste0(output,"_UMAP_Upfemalemodule.pdf"))
FeaturePlot(mySeurat, features="Upfemale.score1", cols = c("light blue", "red"))
dev.off()
