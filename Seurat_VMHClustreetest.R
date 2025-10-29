#!/usr/bin/env Rscript
output <- "Seurat/VMH_SCTpostmerge/VMH_SCTpostmerge_multiclust"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
library(matrixStats)
library(qlcMatrix)
library(igraph)
library(RANN)
library(clustree)
MaleVMH.data <- Read10X(data.dir = "FemalePOAPrimed_Putative2/outs/filtered_feature_bc_matrix")
MaleVMH <- CreateSeuratObject(counts = MaleVMH.data, project = "MaleVMH", min.cells=3, min.features=200)
PrimedVMH.data <- Read10X(data.dir= "FemaleVMHPrimed/outs/filtered_feature_bc_matrix")
PrimedVMH <- CreateSeuratObject(counts = PrimedVMH.data, project = "PrimedVMH", min.cells=3, min.features=200)
MaleVMH$sex <- "Male"
PrimedVMH$sex <- "Female"
mySeurat <- merge(MaleVMH, y=c(PrimedVMH), add.cell.ids=c("MaleVMH","PrimedVMH"), project="Merged")
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress="orig.ident")
mySeurat <- RunPCA(mySeurat)
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
head(mySeurat[[]])
mySeurat <- FindClusters(mySeurat, resolution=c(0.2,0.4,0.6,0.8,1.0,1.2,1.5,2.0), plot.SNN=TRUE, save.SNN=TRUE)
saveRDS(mySeurat, file=(paste0(output, ".rds")))
