#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")
coding <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/vM21_Proteincoding_Genesybols.txt")
coding <- unlist(coding)
coding <- unique(coding)
mCTgenes = data.frame()
mCTcoding = data.frame()
for (i in 0:34)
try({

ident <- i
mCT <- subset(mySeurat, idents=c(i))

mat <- mCT@assays[["RNA"]]@counts

df <- as.data.frame(mat)
dim(df)
collapsed <- rowSums(df)
head(collapsed)
dim(collapsed)
length(collapsed[collapsed>0])
genes <- collapsed[collapsed >0]
length(genes)
genes.df <- as.data.frame(genes)
dim(genes.df)
head(genes.df)
genenames <- rownames(genes.df)
head(genenames)
genenames <- noquote(genenames)
head(genenames)
final <- intersect(genenames,coding)
length(final)
genecount <- length(genenames)
codecount <- length(final)



genefinal <- data.frame(i, genecount)
mCTgenes <- rbind(mCTgenes, genefinal)
codefinal <- data.frame(i, codecount)
mCTcoding = rbind(mCTcoding, codefinal)
})
write.csv(mCTgenes, file=paste0(output,"GenesperCT.csv"))
write.csv(mCTcoding, file=paste0(output,"CodingperCT.csv"))