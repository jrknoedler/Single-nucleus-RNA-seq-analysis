#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/Excitvinhib"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


dsSeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")

excit <- subset(dsSeurat, idents=c("1","3","4","5","6","10","11","12","13","14","17","19","20","21","23","24","25","26","27","29","31","32","33","34","35"), invert=TRUE)
inhib <- subset(dsSeurat, idents=c("1","3","4","5","6","10","11","12","13","14","17","19","20","21","23","24","25","26","27","29","31","32","33","34","35"))
genelist <- read.table("Genelists/POA_MvP_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
data <- dsSeurat@assays[["RNA"]]@counts
data <- as.matrix(data)
excitdata <- excit@assays[["RNA"]]@counts
inhibdata <- inhib@assays[["RNA"]]@counts


excitdata.df <- as.data.frame(t(excitdata))
dim(excitdata.df)
subgroup <- rep(1:3, ceiling(nrow(excitdata.df)/3))
subgroup
dim(subgroup)
excitdata.df$pseudorep <- 0 
size <- nrow(excitdata.df)
size
excitdata.df[,"pseudorep"] <- sample(subgroup, size=size, replace=FALSE)
Pseudorep1 <- excitdata.df[excitdata.df$pseudorep %in% 1,]
dim(Pseudorep1)
Pseudorep2 <- excitdata.df[excitdata.df$pseudorep %in% 2,]
dim(Pseudorep2)
Pseudorep3 <- excitdata.df[excitdata.df$pseudorep %in% 3,]
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
ExcitPseudoreps <- cbind(Pseudobulk1, Pseudobulk2, Pseudobulk3)
head(ExcitPseudoreps)

inhibdata.df <- as.data.frame(t(inhibdata))
dim(inhibdata.df)
subgroup <- rep(1:3, ceiling(nrow(inhibdata.df)/3))
subgroup
dim(subgroup)
inhibdata.df$pseudorep <- 0 
size <- nrow(inhibdata.df)
size
inhibdata.df[,"pseudorep"] <- sample(subgroup, size=size, replace=FALSE)
Pseudorep1 <- inhibdata.df[inhibdata.df$pseudorep %in% 1,]
dim(Pseudorep1)
Pseudorep2 <- inhibdata.df[inhibdata.df$pseudorep %in% 2,]
dim(Pseudorep2)
Pseudorep3 <- inhibdata.df[inhibdata.df$pseudorep %in% 3,]
dim(Pseudorep3)

Pseudorep1 <- t(Pseudorep1)
Pseudorep2 <- t(Pseudorep2)
Pseudorep3 <- t(Pseudorep3)

Pseudobulk1i <- rowSums(Pseudorep1)
Pseudobulk2i <- rowSums(Pseudorep2)
Pseudobulk3i <- rowSums(Pseudorep3)
dim(Pseudobulk1)
dim(Pseudobulk2)
dim(Pseudobulk3)
InhibPseudoreps <- cbind(Pseudobulk1i, Pseudobulk2i, Pseudobulk3i)
head(InhibPseudoreps)

Combined1 <- cbind(ExcitPseudoreps, InhibPseudoreps)
Combined <- as.matrix(Combined1)
Combined <- round(Combined)
cond1 <- "Excitatory"
cond2 <- "Inhibitory"
sampleTable <- data.frame(condition=factor(c(rep(cond1, 3), rep(cond2, 3))))
rownames(sampleTable) <- colnames(Combined)
dds=DESeqDataSetFromMatrix(Combined, sampleTable, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

restargeted <- res[filtered.genelist,]
restargeted$padj <- p.adjust(restargeted$pvalue, method="BH")
write.csv(restargeted, file=paste0(output,"_MvP_pseudobulksDEGpadjust.csv"))
resAll <- results(dds)
table(resAll$padj<0.05)
## Order by adjusted p-value
resAll <- resAll[order(resAll$padj), ]
## Merge with normalized count data
resAlldata <- merge(as.data.frame(resAll), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resAlldata)[1] <- "Gene"
## Write results
write.csv(resAlldata, file=paste0(output, "_Allgenes-ExpDiff.csv"))
