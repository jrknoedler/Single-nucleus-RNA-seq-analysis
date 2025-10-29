#!/usr/bin/env Rscript

output <- "Seurat/VMH_IndependentAnalysis/VMH_Regulons"

library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tidyverse)

mySeurat <- readRDS("/scratch/users/knoedler/Seurat/VMH_IndependentAnalysis/VMHMerged_arcfinalfilter_sexclude_malat1regress_res1.2FINAL.rds")

regulons <- read.csv("/scratch/users/knoedler/pySCENIC_Singularity/VMH_DefaultNES_Final/auc_mtx_filtered.csv", header=TRUE, row.names=1)
dim(regulons)

#regulonlist <- read.csv("/scratch/users/knoedler/pySCENIC_Singularity/BNSTMerged_lowNES/regulonlist2.csv")
#regulonlist <- data.frame(regulonlist)
#regulonlist <- t(regulonlist)
#regulonlist <- unlist(regulonlist)
#regulonlist <- as.matrix(regulonlist)
#regulonlist
idents <- mySeurat@active.ident
idents <- data.frame(idents)
merged <- cbind(regulons, idents)
head(mySeurat[[]])
dim(merged)
#exp <- merged[c(2:328)]
#clusters <- merged[329]

#anno_colors = list(clusters=c("lightgreen","plum","violetred1","lightblue"))

#pdf(file="BNSTRegtest.pdf")
#pheatmap(exp, annotation_row=clusters, cluster_rows=FALSE, cluster_cols=FALSE, color=viridis(500),show_rownames=FALSE, border_color=NA)
#dev.off()

avg <- merged %>% group_by(idents) %>% summarize_all(.funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames(var="idents")
avg <- data.frame(avg)
avg
idents <- avg[1]
#avg %>% remove_rownames %>% column_to_rownames (var="idents")
exp <- avg[c(2:325)]
#exp <- as.matrix(exp)
exp <- scale((exp))
idents
row.names(exp) <- c(0:27)
exp <- data.frame(exp)
pdf(file=paste0(output,"_regbyclust.pdf"), width=40, height=12)
pheatmap(exp, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()
topregs <- colnames(exp)[max.col(exp, ties.method="first")]
topregs
topregs <- unlist(topregs)
topregs <- as.matrix(topregs)
topplot <- exp[,topregs]
pdf(file=paste0(output,"VMH_topregsbyclust.pdf"), width=30, height=12)
pheatmap(topplot, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()
RegDimorphic <- read.table("RegulonsCompiled/VMH/VMHDegulonsSeurat.txt", header=FALSE)
RegDimorphic <- unlist(RegDimorphic)
RegDimorphic <- as.matrix(RegDimorphic)
RegDimorphic
expdimorphic <- exp[,RegDimorphic]

Clustdimorphic <- c("7","22","18","25")
Clustdimorphic <- unlist(Clustdimorphic)
Clustdimorphic <- as.matrix(Clustdimorphic)
dimorphicclusts <- exp[Clustdimorphic,]
dimorphicclusts <- scale(dimorphicclusts)
pdf(file=paste0(output,("VMH_regbyclustdimorphic.pdf", width=20, height=12)
pheatmap(expdimorphic, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()
pdf(file="VMH_dimorphicclusts.pdf", width=40, height=8)
pheatmap(dimorphicclusts, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()
pared <- dimorphicclusts[,RegDimorphic]
pared <- scale(pared)
pdf(file="VMH_paredregulons.pdf")
pheatmap(pared, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()

mySeurat <- AddMetaData(mySeurat, regulons)

pdf(file="VMH_Foxo1Reg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Foxo1..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Gadd45aReg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Gadd45a..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Esr1Reg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Esr1..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_PgrReg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Pgr..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Stat5aReg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Stat5a..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Tead1Reg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Tead1..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_ArReg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Ar..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Bcl6Reg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Bcl6..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_NfibReg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Nfib..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Setbp1Reg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Setbp1..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Sox5Reg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Sox5..."), pt.size=0, ncol=1)
dev.off()
pdf(file="VMH_Nr5a2Reg.pdf")
VlnPlot(mySeurat, idents=c("7","18","22","25"), cols=c("red","blue","red","dark green"), features=c("Nr5a2..."), pt.size=0, ncol=1)
dev.off()