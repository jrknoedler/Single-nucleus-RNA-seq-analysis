#!/usr/bin/env Rscript

sampleID <- "UnprimedBNST"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "ARACHNE/MacoskoLiger/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Femaledata <- Read10X(data.dir = "SCP372/other/BNST_only_Female", gene.column=1)
Female <- CreateSeuratObject(counts = Femaledata, project = "Female", min.cells=3, min.features=200)
Maledata <- Read10X(data.dir = "SCP372/other/BNST_only_Male", gene.column=1)
Male <- CreateSeuratObject(counts = Maledata, project = "Male", min.cells=3, min.features=200)
mySeurat <- merge(Male, y=c(Female), add.cell.ids=c("Male","Female"), project="Merged")
mySeurat.downsample = subset(mySeurat, cells = sample(Cells(mySeurat), 20000))
#mySeurat.downsample <- SCTransform(mySeurat.downsample, verbose=TRUE)


#saveRDS(mySeurat, file=(paste0(output, ".rds")))
mat <- mySeurat.downsample[["RNA"]]@counts

write.csv(mat, file=paste0(output, "expressionmtx.csv"))
