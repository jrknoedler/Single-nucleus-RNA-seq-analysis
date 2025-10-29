#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/Paper_DraftAnalysis/RegressDegome_LessFilt/metaDEGs/PvU_AllDEGs_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")

mySeurat <- subset(mySeurat, idents=c("0","3","8","12","13","14","16","17","20","21","22","28","30","31","32","33","25","4","9"), invert=TRUE)

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)


datas <- SampleUMI(data, max.umi=40000, upsample=FALSE, verbose=TRUE)
datar <- round(datas)

dsSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

hormone <- mySeurat$Hormone
dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
dsSeurat <- subset(dsSeurat, nCount_RNA > 7500)
genelist <- read.table("Genelists/MeA_PvU_1.5cutoff.txt", header=FALSE)
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


Maledata.df <- as.data.frame(t(Maledata))
dim(Maledata.df)
subgroup <- rep(1:3, ceiling(nrow(Maledata.df)/3))
subgroup
dim(subgroup)
Maledata.df$pseudorep <- 0 
size <- nrow(Maledata.df)
size
Maledata.df[,"pseudorep"] <- sample(subgroup, size=size, replace=FALSE)
Pseudorep1 <- Maledata.df[Maledata.df$pseudorep %in% 1,]
dim(Pseudorep1)
Pseudorep2 <- Maledata.df[Maledata.df$pseudorep %in% 2,]
dim(Pseudorep2)
Pseudorep3 <- Maledata.df[Maledata.df$pseudorep %in% 3,]
dim(Pseudorep3)

Pseudorep1 <- t(Pseudorep1)
Pseudorep2 <- t(Pseudorep2)
Pseudorep3 <- t(Pseudorep3)

Pseudobulk1 <- rowSums(Pseudorep1)
Pseudobulk2 <- rowSums(Pseudorep2)
Pseudobulk3 <- rowSums(Pseudorep3)
dim(Pseudobulk1)
dim(Pseudobulk2)
dim(Pseudobulk3)
MalePseudoreps <- cbind(Pseudobulk1, Pseudobulk2, Pseudobulk3)
head(MalePseudoreps)




Femaledata <- Female@assays[["RNA"]]@counts


Femaledata.df <- as.data.frame(t(Femaledata))
dim(Femaledata.df)
subgroupf <- rep(1:3, ceiling(nrow(Femaledata.df)/3))
subgroupf
dim(subgroupf)
Femaledata.df$pseudorep <- 0 
size <- nrow(Femaledata.df)
size
Femaledata.df[,"pseudorep"] <- sample(subgroupf, size=size, replace=FALSE)
Pseudorep1f <- Femaledata.df[Femaledata.df$pseudorep %in% 1,]
Pseudorep2f <- Femaledata.df[Femaledata.df$pseudorep %in% 2,]
Pseudorep3f <- Femaledata.df[Femaledata.df$pseudorep %in% 3,]

Pseudorep1f <- t(Pseudorep1f)
Pseudorep2f <- t(Pseudorep2f)
Pseudorep3f <- t(Pseudorep3f)

Pseudobulk1f <- rowSums(Pseudorep1f)
Pseudobulk2f <- rowSums(Pseudorep2f)
Pseudobulk3f <- rowSums(Pseudorep3f)
dim(Pseudobulk1f)
dim(Pseudobulk2f)
dim(Pseudobulk3f)
FemalePseudoreps <- cbind(Pseudobulk1f, Pseudobulk2f, Pseudobulk3f)
head(FemalePseudoreps)

Combined1 <- cbind(MalePseudoreps, FemalePseudoreps)
Combined <- as.matrix(Combined1)

cond1 <- "Intact"
cond2 <- "Primed"
sampleTable <- data.frame(condition=factor(c(rep(cond1, 3), rep(cond2, 3))))
rownames(sampleTable) <- colnames(Combined)
dds=DESeqDataSetFromMatrix(Combined, sampleTable, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

restargeted <- res[filtered.genelist,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
write.csv(restargeted, file=paste0(output,"_inhib_pseudobulksDEGpadjust.csv"))
