#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/MvP30pcs/BNSTMvP_Cluster3_1.5DESeq2test"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")
mySeurat <- subset(mySeurat, nCount_RNA > 10000)

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)



datas <- SampleUMI(data, max.umi=40000, upsample=FALSE, verbose=TRUE)
datar <- round(datas)

dsSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone


dsSeurat <- AddMetaData(dsSeurat, metadata=idents, col.name="Ident")
dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
Idents(dsSeurat) <- "Ident"
head(dsSeurat[[]])

genelist <- read.table("Genelists/BNST_MvP_1.5cutoff.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
Male <- subset(x=dsSeurat, subset=Hormone=="Intact")
Male <- subset(x=Male, idents=c("3"))
Female <- subset(x=dsSeurat, subset=Hormone=="Primed")
Female <- subset(x=Female, idents=c("3"))
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
write.csv(restargeted, file=paste0(output,"MvP_sDEGpadjust.csv"))

