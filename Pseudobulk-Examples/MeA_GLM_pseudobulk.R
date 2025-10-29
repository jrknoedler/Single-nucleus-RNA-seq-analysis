#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "X:/Seurat/MeA_IndependentAnalysis/MeA_MetaDEGslessfilt"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)
memory.limit(64000)

mySeurat <- readRDS("X:/Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")



data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)


datas <- SampleUMI(data, max.umi=40000, upsample=FALSE, verbose=TRUE)
datar <- round(datas)

dsSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

hormone <- mySeurat$Hormone
dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
dsSeurat <- subset(dsSeurat, nCount_RNA > 5000)
genelist <- read.table("X:/Genelists/MeA_MvP_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.MvP <- intersect(genelist, genes.10x)
filtered.MvP

genelist <- read.table("X:/Genelists/MeA_MvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.MvU <- intersect(genelist, genes.10x)
filtered.MvU

genelist <- read.table("X:/Genelists/MeA_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.PvU <- intersect(genelist, genes.10x)
filtered.PvU

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)

Male <- subset(x=dsSeurat, subset=Hormone=="Intact")
head(Male[[]])
Female <- subset(x=dsSeurat, subset=Hormone=="Primed")
Diestrus <- subset(x=dsSeurat, subset=Hormone=="Unprimed")
Maledata <- Male@assays[["RNA"]]@counts


Maledata.df <- as.data.frame(t(Maledata))
dim(Maledata.df)
subgroup <- rep(1:5, ceiling(nrow(Maledata.df)/5))
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
Pseudorep4 <- Maledata.df[Maledata.df$pseudorep %in% 4,]
Pseudorep5 <- Maledata.df[Maledata.df$pseudorep %in% 5,]

Pseudorep1 <- t(Pseudorep1)
Pseudorep2 <- t(Pseudorep2)
Pseudorep3 <- t(Pseudorep3)
Pseudorep4 <- t(Pseudorep4)
Pseudorep5 <- t(Pseudorep5)


Pseudobulk1 <- rowSums(Pseudorep1)
Pseudobulk2 <- rowSums(Pseudorep2)
Pseudobulk3 <- rowSums(Pseudorep3)
Pseudobulk4 <- rowSums(Pseudorep4)
Pseudobulk5 <- rowSums(Pseudorep5)

dim(Pseudobulk1)
dim(Pseudobulk2)
dim(Pseudobulk3)
MalePseudoreps <- cbind(Pseudobulk1, Pseudobulk2, Pseudobulk3, Pseudobulk4, Pseudobulk5)
head(MalePseudoreps)




Femaledata <- Female@assays[["RNA"]]@counts


Femaledata.df <- as.data.frame(t(Femaledata))
dim(Femaledata.df)
subgroupf <- rep(1:5, ceiling(nrow(Femaledata.df)/5))
subgroupf
dim(subgroupf)
Femaledata.df$pseudorep <- 0 
size <- nrow(Femaledata.df)
size
Femaledata.df[,"pseudorep"] <- sample(subgroupf, size=size, replace=FALSE)
Pseudorep1f <- Femaledata.df[Femaledata.df$pseudorep %in% 1,]
Pseudorep2f <- Femaledata.df[Femaledata.df$pseudorep %in% 2,]
Pseudorep3f <- Femaledata.df[Femaledata.df$pseudorep %in% 3,]
Pseudorep4f <- Femaledata.df[Femaledata.df$pseudorep %in% 4,]
Pseudorep5f <- Femaledata.df[Femaledata.df$pseudorep %in% 5,]

Pseudorep1f <- t(Pseudorep1f)
Pseudorep2f <- t(Pseudorep2f)
Pseudorep3f <- t(Pseudorep3f)
Pseudorep4f <- t(Pseudorep4f)
Pseudorep5f <- t(Pseudorep5f)

Pseudobulk1f <- rowSums(Pseudorep1f)
Pseudobulk2f <- rowSums(Pseudorep2f)
Pseudobulk3f <- rowSums(Pseudorep3f)
Pseudobulk4f <- rowSums(Pseudorep4f)
Pseudobulk5f <- rowSums(Pseudorep5f)

dim(Pseudobulk1f)
dim(Pseudobulk2f)
dim(Pseudobulk3f)
FemalePseudoreps <- cbind(Pseudobulk1f, Pseudobulk2f, Pseudobulk3f, Pseudobulk4f, Pseudobulk5f)
head(FemalePseudoreps)


Unprimeddata <- Diestrus@assays[["RNA"]]@counts

Unprimeddata.df <- as.data.frame(t(Unprimeddata))

subgroupu <- rep(1:5, ceiling(nrow(Unprimeddata.df)/5))

dim(subgroupu)
Unprimeddata.df$pseudorep <- 0 
size <- nrow(Unprimeddata.df)
size
Unprimeddata.df[,"pseudorep"] <- sample(subgroupu, size=size, replace=FALSE)
Pseudorep1u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 1,]
Pseudorep2u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 2,]
Pseudorep3u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 3,]
Pseudorep4u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 4,]
Pseudorep5u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 5,]

Pseudorep1u <- t(Pseudorep1u)
Pseudorep2u <- t(Pseudorep2u)
Pseudorep3u <- t(Pseudorep3u)
Pseudorep4u <- t(Pseudorep4u)
Pseudorep5u <- t(Pseudorep5u)

Pseudobulk1u <- rowSums(Pseudorep1u)
Pseudobulk2u <- rowSums(Pseudorep2u)
Pseudobulk3u <- rowSums(Pseudorep3u)
Pseudobulk4u <- rowSums(Pseudorep4u)
Pseudobulk5u <- rowSums(Pseudorep5u)

UnprimedPseudoreps <- cbind(Pseudobulk1u, Pseudobulk2u, Pseudobulk3u, Pseudobulk4u, Pseudobulk5u)




Combined1 <- cbind(MalePseudoreps, FemalePseudoreps, UnprimedPseudoreps)
Combined <- as.matrix(Combined1)



cond1 <- "Intact"
cond2 <- "Primed"
cond3 <- "Unprimed"

sampleTable <- data.frame(condition=factor(c(rep(cond1, 5), rep(cond2, 5), rep(cond3,5))))
rownames(sampleTable) <- colnames(Combined)
dds=DESeqDataSetFromMatrix(Combined, sampleTable, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

resMalevPrimed <- results(dds, contrast=c("condition",cond1,cond2))
table(resMalevPrimed$padj<0.05)
## Order by adjusted p-value
restargeted <- resMalevPrimed[filtered.MvP,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
table(restargeted$padj<0.05)
write.csv(restargeted, file=paste0(output,"MvP_GLMpseudo.csv"))

resMalevUnprimed <- results(dds, contrast=c("condition",cond1,cond3))
table(resMalevUnprimed$padj<0.05)
## Order by adjusted p-value
restargeted <- resMalevUnprimed[filtered.MvU,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
table(restargeted$padj<0.05)
write.csv(restargeted, file=paste0(output,"MvU_GLMpseudo.csv"))

resPrimedvUnprimed <- results(dds, contrast=c("condition",cond2,cond3))
table(resPrimedvUnprimed$padj<0.05)
## Order by adjusted p-value
restargeted <- resPrimedvUnprimed[filtered.PvU,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
table(restargeted$padj<0.05)
write.csv(restargeted, file=paste0(output,"PvU_GLMpseudo.csv"))


