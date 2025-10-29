#!/usr/bin/env Rscript
output <- "Seurat/VMH_IndependentAnalysis/PrimedSDEG_20PCS_Plottest"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1_20pcs_Primed.rds")
head(mySeurat[[]])
mySeurat
genelist <- read.table("Genelists/VMH_MalevPrimedAll.txt")
unlist(genelist)
genelist <- as.matrix(genelist)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
SFARIgenes <- read.table("Genelists/SFARI.txt")
unlist(SFARIgenes)
SFARIgenes <- as.matrix(SFARIgenes)
filtered.SFARI <- intersect(SFARIgenes, genes.10x) 
Upmale <- read.table("Genelists/VMH_MalevPrimed_Maleup.txt")
unlist(Upmale)
Upmale <- as.matrix(Upmale)
filtered.Upmale <- intersect(Upmale, genes.10x)
Upfemale <- read.table("Genelists/VMH_MalevPrimed_Primedup.txt")
unlist(Upfemale)
Upfemale <- as.matrix(Upfemale)
filtered.Upfemale <- intersect(Upfemale, genes.10x)
UpPvU <- read.table("Genelists/VMH_PvU_UpPrimed.txt")
unlist(UpPvU)
UpPvU <- as.matrix(UpPvU)
filtered.UpPvU <- intersect(UpPvU, genes.10x)
UpUvP <- read.table("Genelists/VMH_PvU_UpUnprimed.txt")
unlist(UpUvP)
UpUvP <- as.matrix(UpUvP)
filtered.UpUvP <- intersect(UpUvP, genes.10x)
UpMaleCombined <- read.table("Genelists/VMH_UpMaleCombined.txt")
unlist(UpMaleCombined)
UpMaleCombined <- as.matrix(UpMaleCombined)
UpMaleCombined.filtered <- intersect(UpMaleCombined, genes.10x)
UpPrimedCombined <- read.table("Genelists/VMH_UpPrimedCombined.txt")
unlist(UpPrimedCombined)
UpPrimedCombined <- as.matrix(UpPrimedCombined)
UpPrimedCombined.filtered <- intersect(UpPrimedCombined, genes.10x)
UpUnprimedCombined <- read.table("Genelists/VMH_UpUnprimedCombined.txt")
unlist(UpUnprimedCombined)
UpUnprimedCombined <- as.matrix(UpUnprimedCombined)
UpUnprimedCombined.filtered <- intersect(UpUnprimedCombined, genes.10x)
mySeurat[["percent.SDEG"]] <- PercentageFeatureSet(object=mySeurat, features=filtered.genelist)
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output,"_UMAP_pct_SDEG.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="percent.SDEG", cols = c("light blue", "red"))
dev.off()
pdf(paste0(output,"_UMAP_Esr1.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Esr1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Cyp19a1.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Cyp19a1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Ar.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Ar", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Cckar.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Cckar", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Phf21b.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Phf21b", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Setdb2.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Setdb2", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Pdyn.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Pdyn", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Gldn.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Gldn", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Mical2.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Mical2", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Tnxb.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Tnxb", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Ezr.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Ezr", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Maml2.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Maml2", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Iqschfp.pdf"))
FeaturePlot(mySeurat, reduction="umap", features="Iqschfp", cols = c("light blue", "red"), order=TRUE)
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
mySeurat <- AddModuleScore(
object = mySeurat,
features = list(filtered.UpPvU),
name = 'UpPvU.score'
)
mySeurat <- AddModuleScore(
object = mySeurat,
features = list(filtered.UpPvU),
name = 'UpUvP.score'
)
mySeurat <- AddModuleScore(
object = mySeurat,
features=list(UpMaleCombined.filtered),
name = 'Upmalecombined.score'
)
mySeurat <- AddModuleScore(
object = mySeurat,
features=list(UpPrimedCombined.filtered),
name = 'UpPrimedcombined.score'
)
mySeurat <- AddModuleScore(
object = mySeurat,
features=list(UpUnprimedCombined.filtered),
name = 'UpUnprimedcombined.score'
)
head(mySeurat[[]])
pdf(paste0(output,"_UMAP_SDEGmodule.pdf"))
FeaturePlot(mySeurat, features="SDEG.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_SFARImodule.pdf"))
FeaturePlot(mySeurat, features="SFARI.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Upmalemodule.pdf"))
FeaturePlot(mySeurat, features="Upmale.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_UpPvUmodule.pdf"))
FeaturePlot(mySeurat, features="UpPvU.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_UpUvPmodule.pdf"))
FeaturePlot(mySeurat, features="UpUvP.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_UpMaleCombined.pdf"))
FeaturePlot(mySeurat, features="Upmalecombined.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_UpPrimedCombined.pdf"))
FeaturePlot(mySeurat, features="UpPrimedcombined.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_UpUnprimedCombined.pdf"))
FeaturePlot(mySeurat, features="UpUnprimedcombined.score1", cols = c("light blue", "red"), order=TRUE)
dev.off()
