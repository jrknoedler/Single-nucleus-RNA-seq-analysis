#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MalePOA_Preprocessing/PrimedPOA_preprocessing_test"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Seurat.data <- Read10X(data.dir = "FemalePOAPrimed/outs/filtered_feature_bc_matrix")
mySeurat <- CreateSeuratObject(counts = Seurat.data, project = sampleID, min.cells=3, min.features=200)
antisense <- read.table("Genelists/antisensenames.txt", header=FALSE)
antisense <- unlist(antisense)
head(antisense)
genes <- rownames(mySeurat)
genes <- as.matrix(genes)
head(genes)
antisense.present <- intersect(antisense, genes)
antisense.present
mySeurat[["percent.antisense"]] <- PercentageFeatureSet(mySeurat, features=antisense.present)
pdf(file=paste0(output, "_antisensecounts.pdf"))
VlnPlot(mySeurat, features = c("percent.antisense"))
dev.off()
