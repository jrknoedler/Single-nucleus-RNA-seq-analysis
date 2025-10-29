#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/GLM_MetaDEGs10repBY"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)

memory.limit(64000)
mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")



data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)


datas <- SampleUMI(data, max.umi=40000, upsample=FALSE, verbose=TRUE)
datar <- round(datas)

dsSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

hormone <- mySeurat$Hormone
dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
dsSeurat <- subset(dsSeurat, nCount_RNA > 5000)
genelist <- read.table("Genelists/BNST_MvP_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- unique(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.MvP <- intersect(genelist, genes.10x)


genelist <- read.table("Genelists/BNST_MvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.MvU <- intersect(genelist, genes.10x)
filtered.MvU

genelist <- read.table("Genelists/BNST_PvU_1.5cutoff.txt", header=FALSE)
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
subgroup <- rep(1:10, ceiling(nrow(Maledata.df)/10))
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
Pseudorep6 <- Maledata.df[Maledata.df$pseudorep %in% 6,]
Pseudorep7 <- Maledata.df[Maledata.df$pseudorep %in% 7,]
Pseudorep8 <- Maledata.df[Maledata.df$pseudorep %in% 8,]
Pseudorep9 <- Maledata.df[Maledata.df$pseudorep %in% 9,]
Pseudorep10 <- Maledata.df[Maledata.df$pseudorep %in% 10,]



Pseudorep1 <- t(Pseudorep1)
Pseudorep2 <- t(Pseudorep2)
Pseudorep3 <- t(Pseudorep3)
Pseudorep4 <- t(Pseudorep4)
Pseudorep5 <- t(Pseudorep5)
Pseudorep6 <- t(Pseudorep6)
Pseudorep7 <- t(Pseudorep7)
Pseudorep8 <- t(Pseudorep8)
Pseudorep9 <- t(Pseudorep9)
Pseudorep10 <- t(Pseudorep10)



Pseudobulk1 <- rowSums(Pseudorep1)
Pseudobulk2 <- rowSums(Pseudorep2)
Pseudobulk3 <- rowSums(Pseudorep3)
Pseudobulk4 <- rowSums(Pseudorep4)
Pseudobulk5 <- rowSums(Pseudorep5)
Pseudobulk6 <- rowSums(Pseudorep6)
Pseudobulk7 <- rowSums(Pseudorep7)
Pseudobulk8 <- rowSums(Pseudorep8)
Pseudobulk9 <- rowSums(Pseudorep9)
Pseudobulk10 <- rowSums(Pseudorep10)


dim(Pseudobulk1)
dim(Pseudobulk2)
dim(Pseudobulk3)
MalePseudoreps <- cbind(Pseudobulk1, Pseudobulk2, Pseudobulk3, Pseudobulk4, Pseudobulk5, Pseudobulk6, Pseudobulk7, Pseudobulk8, Pseudobulk9, Pseudobulk10)
head(MalePseudoreps)




Femaledata <- Female@assays[["RNA"]]@counts


Femaledata.df <- as.data.frame(t(Femaledata))
dim(Femaledata.df)
subgroupf <- rep(1:10, ceiling(nrow(Femaledata.df)/10))
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
Pseudorep6f <- Femaledata.df[Femaledata.df$pseudorep %in% 6,]
Pseudorep7f <- Femaledata.df[Femaledata.df$pseudorep %in% 7,]
Pseudorep8f <- Femaledata.df[Femaledata.df$pseudorep %in% 8,]
Pseudorep9f <- Femaledata.df[Femaledata.df$pseudorep %in% 9,]
Pseudorep10f <- Femaledata.df[Femaledata.df$pseudorep %in% 10,]


Pseudorep1f <- t(Pseudorep1f)
Pseudorep2f <- t(Pseudorep2f)
Pseudorep3f <- t(Pseudorep3f)
Pseudorep4f <- t(Pseudorep4f)
Pseudorep5f <- t(Pseudorep5f)
Pseudorep6f <- t(Pseudorep6f)
Pseudorep7f <- t(Pseudorep7f)
Pseudorep8f <- t(Pseudorep8f)
Pseudorep9f <- t(Pseudorep9f)
Pseudorep10f <- t(Pseudorep10f)


Pseudobulk1f <- rowSums(Pseudorep1f)
Pseudobulk2f <- rowSums(Pseudorep2f)
Pseudobulk3f <- rowSums(Pseudorep3f)
Pseudobulk4f <- rowSums(Pseudorep4f)
Pseudobulk5f <- rowSums(Pseudorep5f)
Pseudobulk6f <- rowSums(Pseudorep6f)
Pseudobulk7f <- rowSums(Pseudorep7f)
Pseudobulk8f <- rowSums(Pseudorep8f)
Pseudobulk9f <- rowSums(Pseudorep9f)
Pseudobulk10f <- rowSums(Pseudorep10f)



dim(Pseudobulk1f)
dim(Pseudobulk2f)
dim(Pseudobulk3f)
FemalePseudoreps <- cbind(Pseudobulk1f, Pseudobulk2f, Pseudobulk3f, Pseudobulk4f, Pseudobulk5f, Pseudobulk6f, Pseudobulk7f, Pseudobulk8f, Pseudobulk9f, Pseudobulk10f)
head(FemalePseudoreps)


Unprimeddata <- Diestrus@assays[["RNA"]]@counts

Unprimeddata.df <- as.data.frame(t(Unprimeddata))

subgroupu <- rep(1:10, ceiling(nrow(Unprimeddata.df)/10))

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
Pseudorep6u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 6,]
Pseudorep7u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 7,]
Pseudorep8u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 8,]
Pseudorep9u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 9,]
Pseudorep10u <- Unprimeddata.df[Unprimeddata.df$pseudorep %in% 10,]




Pseudorep1u <- t(Pseudorep1u)
Pseudorep2u <- t(Pseudorep2u)
Pseudorep3u <- t(Pseudorep3u)
Pseudorep4u <- t(Pseudorep4u)
Pseudorep5u <- t(Pseudorep5u)
Pseudorep6u <- t(Pseudorep6u)
Pseudorep7u <- t(Pseudorep7u)
Pseudorep8u <- t(Pseudorep8u)
Pseudorep9u <- t(Pseudorep9u)
Pseudorep10u <- t(Pseudorep10u)


Pseudobulk1u <- rowSums(Pseudorep1u)
Pseudobulk2u <- rowSums(Pseudorep2u)
Pseudobulk3u <- rowSums(Pseudorep3u)
Pseudobulk4u <- rowSums(Pseudorep4u)
Pseudobulk5u <- rowSums(Pseudorep5u)
Pseudobulk6u <- rowSums(Pseudorep6u)
Pseudobulk7u <- rowSums(Pseudorep7u)
Pseudobulk8u <- rowSums(Pseudorep8u)
Pseudobulk9u <- rowSums(Pseudorep9u)
Pseudobulk10u <- rowSums(Pseudorep10u)


UnprimedPseudoreps <- cbind(Pseudobulk1u, Pseudobulk2u, Pseudobulk3u, Pseudobulk4u, Pseudobulk5u, Pseudobulk6u, Pseudobulk7u, Pseudobulk8u, Pseudobulk9u, Pseudobulk10u)




Combined1 <- cbind(MalePseudoreps, FemalePseudoreps, UnprimedPseudoreps)
Combined <- as.matrix(Combined1)



cond1 <- "Intact"
cond2 <- "Primed"
cond3 <- "Unprimed"

sampleTable <- data.frame(condition=factor(c(rep(cond1, 10), rep(cond2, 10), rep(cond3,10))))
rownames(sampleTable) <- colnames(Combined)
dds=DESeqDataSetFromMatrix(Combined, sampleTable, design=~condition)


dds <- DESeq(dds)
res <- results(dds)

resMalevPrimed <- results(dds, contrast=c("condition",cond1,cond2))
table(resMalevPrimed$padj<0.1)
## Order by adjusted p-value
restargeted <- resMalevPrimed[filtered.MvP,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
table(restargeted$padj<0.1)
write.csv(restargeted, file=paste0(output,"MvP_GLMpseudo.csv"))

resMalevUnprimed <- results(dds, contrast=c("condition",cond1,cond3))
table(resMalevUnprimed$padj<0.1)
## Order by adjusted p-value
restargeted <- resMalevUnprimed[filtered.MvU,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
table(restargeted$padj<0.05)
write.csv(restargeted, file=paste0(output,"MvU_GLMpseudo.csv"))

resPrimedvUnprimed <- results(dds, contrast=c("condition",cond2,cond3))
table(resPrimedvUnprimed$padj<0.1)
## Order by adjusted p-value
restargeted <- resPrimedvUnprimed[filtered.PvU,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
table(restargeted$padj<0.1)
write.csv(restargeted, file=paste0(output,"PvU_GLMpseudo.csv"))




