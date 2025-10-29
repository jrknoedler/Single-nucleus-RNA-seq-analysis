#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome/MetaDEGs/PvU_AllDEGs_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")

mySeurat <- subset(mySeurat, idents=c("0","2","7","8","9","15","16","18","22","28","30","36","37","38"))

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)


datas <- SampleUMI(data, max.umi=40000, upsample=FALSE, verbose=TRUE)
datar <- round(datas)

dsSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

hormone <- mySeurat$Hormone
dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
dsSeurat <- subset(dsSeurat, nCount_RNA > 7500)
genelist <- read.table("Genelists/POA_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)

Male <- subset(x=dsSeurat, subset=Hormone=="Primed")
head(Male[[]])
Female <- subset(x=dsSeurat, subset=Hormone=="Unprimed")
Maledata <- Male@assays[["RNA"]]@counts
Male.mat <- as.matrix(Maledata)
MaleNum <- ncol(Male.mat)
Femaledata <- Female@assays[["RNA"]]@counts
Female.mat <- as.matrix(Femaledata)
FemaleNum <- ncol(Female.mat)
Combined <- cbind(Male.mat, Female.mat)
ncol(Combined)
cond1 <- "Intact"
cond2 <- "Primed"
sampleTable <- data.frame(condition=factor(c(rep(cond1, MaleNum), rep(cond2, FemaleNum))))
rownames(sampleTable) <- colnames(Combined)
dds=DESeqDataSetFromMatrix(Combined, sampleTable, design=~condition)
dds <- DESeq(dds)
res <- results(dds)
restargeted <- res[filtered.genelist,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
write.csv(restargeted, file=paste0(output,"_excitsDEGpadjust.csv"))
