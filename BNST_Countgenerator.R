#!/usr/bin/env Rscript

library(Seurat)

Primed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_Primed_40pcs_res1.5_lowthresh_Primed.rds")
data1 <- Primed@assays[["RNA"]]@counts
data1 <- as.matrix(data1)
data1 <- t(data1)
write.table(data1, file="Seurat_Countmatrices/BNST/BNST_Primed_counts.tsv",quote=FALSE, sep="\t", col.names=TRUE)
write.table(Primed@active.ident, file="Seurat_Countmatrices/BNST/BNST_Primed_Idents.tsv", quote=FALSE, sep="\t", col.names=TRUE)

Unprimed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered2_Unprimed_40pcs_nothresh_res1.5__Unprimed.rds")
data2 <- Unprimed@assays[["RNA"]]@counts
data2 <- as.matrix(data2)
data2 <- t(data2)
write.table(data2, file="Seurat_Countmatrices/BNST/BNST_Unprimed_counts.tsv",quote=FALSE, sep="\t", col.names=TRUE)
write.table(Unprimed@active.ident, file="Seurat_Countmatrices/BNST/BNST_Unprimed_Idents.tsv", quote=FALSE, sep="\t", col.names=TRUE)