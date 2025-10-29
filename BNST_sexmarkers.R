#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_Plotsupdate_3_9"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(pheatmap)
library(viridis)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Figs5heatmapsfinal_"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/BNST_MvP_1.5cutoff.txt")
eDEG <- read.table("Genelists/BNST_PvU_1.5cutoff.txt")
sDEG2 <- read.table("Genelists/BNST_MvU_1.5cutoff.txt")
unlist(eDEG)
unlist(sDEG2)
unlist(Dimorphic)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered.m <- as.matrix(eDEG.filtered)

sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered.m <- as.matrix(sDEG2.filtered)

Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered.m <- as.matrix(Dimorphic.filtered)


AllsDegs <- read.table("topGO/Total_ByRegion/BNST_Genesonly.txt")
AllsDegs <- unlist(AllsDegs)
AllsDegs.m <- as.matrix(AllsDegs)
AllsDegs.filtered <- intersect(AllsDegs, genes.10x)
AllsDegs.filtered <- as.matrix(AllsDegs.filtered)

mySeurat <- NormalizeData(mySeurat)
MvP.markers <- FindAllMarkers(mySeurat, features=Dimorphic.filtered, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.0.223, test.use="MAST")
head(MvP.markers)
Sig.MvP <- MvP.markers[MvP.markers$p_val_adj < 0.05,]
head(Sig.MvP)
genesMvP <- Sig.MvP$gene
head(genesMvP)
uniqueMvP <- names(table(genesMvP)[table(genesMvP)==1])
uniqueMvP
uniqueMvP <- unlist(uniqueMvP)
uniqueMvP <- as.matrix(uniqueMvP)
uniqueclustsMvP <- Sig.MvP[uniqueMvP,]
uniqueclustsMvP
write.csv(uniqueclustsMvP, file=paste0(output,"MvP_uniquemarkerswithclusters.csv"))
mySeurat <- ScaleData(mySeurat)
Avg <- AverageExpression(mySeurat, assays=c("RNA"), features=uniqueMvP)
Avg
write.csv(Avg, file=paste0(output,"MvPUniquemarkers_AvgExpression.csv"))
names(Avg) = gsub(pattern="RNA.", replacement="", x=names(Avg))
Avg
Avg2 <- Avg$RNA
Avg2
MarkedmCTsMvP <- uniqueclustsMvP$cluster
MarkedmCTsMvP <- unique(MarkedmCTsMvP)
MarkedAvg <- Avg2[,MarkedmCTsMvP]
MarkedAvg
MarkedAvg <- scale(t(MarkedAvg))
pdf(file=paste0(output,"MvPSigmarkerDEGs.pdf"))
pheatmap(MarkedAvg, cluster_rows=FALSE, cluster_cols=FALSE, col=magma(2000))
dev.off()
write.csv(MarkedAvg, file=paste0(output,"SexMarkersscaled.csv"))


