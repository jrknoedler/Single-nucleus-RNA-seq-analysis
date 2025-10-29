#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Paper1stpass/20pcs/BNST_IntactvPrimed_Cluster_DESeq2downsample"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_20pcs.rds")

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)

datar <- round(data)

datars <- SampleUMI(datar, max.umi=20000, upsample=FALSE, verbose=TRUE)

RoundedSeurat <- CreateSeuratObject(datars, project="IntegerMatrix", assay="RNA")

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone

RoundedSeurat <- AddMetaData(RoundedSeurat, metadata=idents, col.name="Ident")
RoundedSeurat <- AddMetaData(RoundedSeurat, metadata=hormone, col.name="Hormone")
Idents(RoundedSeurat) <- "Ident"
head(RoundedSeurat[[]])


RoundedSeurat$celltype.Hormone <- paste(Idents(RoundedSeurat), RoundedSeurat$Hormone, sep="_")
RoundedSeurat$celltype <- Idents(RoundedSeurat)
Idents(RoundedSeurat) <- "celltype.Hormone"
for (i in 0:31){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i,"_Primed")
sex.dimorphism <- FindMarkers(RoundedSeurat, slot="counts", ident.1 = ident1, ident.2=ident2, verbose=TRUE, test.use="DESeq2")
write.csv(sex.dimorphism, file=paste0(output,i,"_dimorphism.csv"))
})
}
