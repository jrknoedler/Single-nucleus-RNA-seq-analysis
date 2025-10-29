#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"


library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Figs5heatmapsfinal_markers_alt"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_sexclude_malat1regress_30pcsres1.2drd3filt3_35pcs.rds")
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
AllsDegs.filtered.m <- as.matrix(AllsDegs.filtered)
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
#Markers <- FindAllMarkers(mySeurat, assay="RNA", only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.20)
#write.csv(Markers, file=paste0(output,"allMarkerswtfwilcox.csv"))
All.markers <- FindAllMarkers(mySeurat, assay="RNA", features=AllsDegs.filtered, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.223, test.use="MAST")
write.csv(All.markers, file=paste0(output,"_RNA.csv"))
#saveRDS(Markers, file=paste0(output,"precomputed.rds"))
Sig.Markers <- All.markers[All.markers$p_val_adj < 0.05,]
Sig.Markers <- Sig.Markers[Sig.Markers$pct.2 < 0.25,]
head(Sig.Markers)
siggenes <- Sig.Markers$gene
head(siggenes)
uniquesig <- names(table(siggenes)[table(siggenes)==1])
uniquesig
uniquesig <- unlist(uniquesig)
uniquesig.m <- as.matrix(uniquesig)
uniqueclustssig <- Sig.Markers[uniquesig.m,]
uniqueclustssig
MvP <- uniqueclustssig[Dimorphic.filtered.m,]
MvU <- uniqueclustssig[sDEG2.filtered.m,]
PvU <- uniqueclustssig[eDEG.filtered.m,]
MarkedMvP <- MvP$cluster
MarkedMvP <- unique(MarkedMvP)
MarkedMvP <- na.omit(MarkedMvP)
MarkedMvU <- MvU$cluster
MarkedMvU <- unique(MarkedMvU)
MarkedPvU <- PvU$cluster
MarkedPvU <- unique(PvU)
write.csv(uniqueclustssig, file=paste0(output,"MvP_uniquemarkerswithclusters.csv"))
Avg <- AverageExpression(mySeurat, assays=c("RNA"), features=uniquesig)
Avg
write.csv(Avg, file=paste0(output,"allsigUniquemarkers_AvgExpression.csv"))
names(Avg) = gsub(pattern="RNA.", replacement="", x=names(Avg))
Avg
Avg2 <- Avg$RNA
Avg2
MarkedmCTssig <- uniqueclustssig$cluster
MarkedmCTssig  <- unique(MarkedmCTssig)
MarkedAvg <- Avg2[,MarkedmCTssig]
MarkedAvg
MarkedAvg <- scale(t(MarkedAvg))
pdf(file=paste0(output,"allSigmarkerDEGs.pdf"))
pheatmap(MarkedAvg, col=magma(2000))
dev.off()
write.csv(MarkedAvg, file=paste0(output,"SexMarkersscaled.csv"))
MvPavg <- Avg2[Dimorphic.filtered.m,]
pdf(file=paste0(output, "Dotplot25.pdf"), width=20)
DotPlot(mySeurat, features=uniquesig, dot.min=0.25)
dev.off()
pdf(file=paste0(output, "Dotplot10.pdf"), width=20)
DotPlot(mySeurat, features=uniquesig, dot.min=0.10)
dev.off()
