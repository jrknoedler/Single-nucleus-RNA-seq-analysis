#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_RegulonsFebredo"

library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(tidyverse)

mySeurat <- readRDS("/scratch/users/knoedler/Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered2_sexclude_malat1regress_30pcsres1.2FINAL.rds")

regulons <- read.csv("/scratch/users/knoedler/pySCENIC_Singularity/BNSTMerged_lowNES_Final/auc_mtx_filtered.csv", header=TRUE, row.names=1)
#mySeurat <- AddMetaData(mySeurat, regulons)
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
#exp <- merged[c(2:410)]
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
exp <- avg[c(2:328)]
#exp <- as.matrix(exp)
exp <- scale((exp))
idents
row.names(exp) <- c(0:34)
exp <- data.frame(exp)
pdf(file="BNST_regbyclust.pdf", width=40, height=12)
pheatmap(exp, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()
topregs <- colnames(exp)[max.col(exp, ties.method="first")]
topregs
topregs <- unlist(topregs)
topregs <- as.matrix(topregs)
topplot <- exp[,topregs]
pdf(file=paste0(output,"_topregsbyclust.pdf"), width=40, height=15)
pheatmap(topplot, show_rownames=T, border_color=NA, color=viridis(2000))
dev.off()
RegDimorphic <- read.table("RegulonsCompiled/BNST/RegDegs.txt", header=FALSE)
RegDimorphic <- unlist(RegDimorphic)
RegDimorphic <- as.matrix(RegDimorphic)
RegDimorphic
expdimorphic <- exp[,RegDimorphic]

pdf(file=paste0(output,"_regbyclustdimorphic.pdf"), width=20, height=12)
pheatmap(expdimorphic, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()

DegulonDimorphic <- read.table("RegulonsCompiled/BNST/10DEGs.txt", header=FALSE)
DegulonDimorphic <- unlist(DegulonDimorphic)
DegulonDimorphic <- as.matrix(DegulonDimorphic)
DegulonDimorphic
expdegdimorphic <- exp[,DegulonDimorphic]

pdf(file="BNST_regbyclustdegulons_targeted.pdf", width=20, height=12)
pheatmap(expdegdimorphic, show_rownames=T, border_color=NA, color=viridis(500))
dev.off()