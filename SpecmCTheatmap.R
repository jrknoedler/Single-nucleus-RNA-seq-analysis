#!/usr/bin/env Rscript

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)
output <- "Heatmaps/FullData_Uniquemarkers"
BNST <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")

MeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")

POA <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered6_redo.rds")

VMH <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")

BNST$celltype.region <- paste(Idents(BNST), "BNST", sep="_")
BNST$celltype <- Idents(BNST)
Idents(BNST) <- "celltype.region"

MeA$celltype.region <- paste(Idents(MeA), "MeA", sep="_")
MeA$celltype <- Idents(MeA)
Idents(MeA) <- "celltype.region"

POA$celltype.region <- paste(Idents(POA), "POA", sep="_")
POA$celltype <- Idents(POA)
Idents(POA) <- "celltype.region"

VMH$celltype.region <- paste(Idents(VMH), "VMH", sep="_")
VMH$celltype <- Idents(VMH)
Idents(VMH) <- "celltype.region"

mySeurat <- merge(BNST, y=c(MeA,POA,VMH), add.cell.ids=c("BNST","MeA","POA","VMH"), project="Merged")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)

DEGs <- read.csv("Genelists/All_uniquemarkers.csv", header=FALSE)
DEGs <- unlist(DEGs)
DEGs <- as.matrix(DEGs)
Avg <- AverageExpression(mySeurat, assay="RNA")
Avg2 <- Avg$RNA
DEGavg <- Avg2[DEGs,]
write.csv(DEGavg, file=paste0(output,"DEGmatrix.csv"))
scaled <- scale(t(DEGavg))
pdf(file=paste0(output,"Rawheatmap.pdf"), width=40, height=40)
pheatmap(scaled, col=magma(2000))
dev.off()
write.csv(scaled, file=paste0(output,"scaledDEGmatrix.csv"))
clusters <- c("9_BNST","17_BNST","19_BNST","20_BNST","32_BNST","33_BNST","35_BNST","32_MeA","5_MeA","4_MeA","10_MeA","14_MeA","30_POA","0_POA","28_POA","1_POA","2_POA","22_POA","27_POA","20_POA","31_POA","12_POA","18_POA","26_POA","19_VMH","24_VMH","20_VMH","23_VMH")
clusters.df <- data.frame(clusters)
scaledsub <- scaled[rownames(clusters) %in% clusters,]
pdf(file=paste0(output,"Subsetheatmap.pdf"), width=40,height=40)
pheatmap(scaledsub, col=magma(2000))
dev.off()
write.csv(scaledsub, file=paste0(output,"subsetscaled.csv"))
DEGmCTavg <- Avg2[DEGs,clusters]
scaled2 <- scale(t(DEGmCTavg))
pdf(file=paste0(output,"Scaledsubset.pdf"), width=40, height=40)
pheatmap(scaled2, cluster_rows=FALSE, col=magma(2000))
dev.off()