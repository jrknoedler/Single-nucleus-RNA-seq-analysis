#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_MalevUnprimedDiff/Scatter/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/BNST_MalevUnprimed_Exploratory/BNST_MaleUnprimed_doubletsremoved_filtered1.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
head(mySeurat[[]])
mySeurat <- NormalizeData(mySeurat)
MaleUP <- read.table("Genelists/BNST_MalevUnprimed_Upmale.txt", header=FALSE)
MaleUP <- unlist(MaleUP)
FemaleUP <- read.table("Genelists/BNST_MalevUnprimed_UpUnprimed.txt", header=FALSE)
FemaleUP <- unlist(FemaleUP)
for (i in 0:42){
try({
theme_set(theme_cowplot())
cluster <- subset(mySeurat, idents = c(i))
Idents(cluster) <- "sex"
genes.10x <- (x=rownames(x=cluster))
filtered.MaleUP <- intersect(MaleUP, genes.10x)
filtered.FemaleUP <- intersect(FemaleUP, genes.10x)
unchanged1 <- setdiff(genes.10x, MaleUP)
unchanged2 <- setdiff(unchanged1, FemaleUP)
avg.exp <- log1p(AverageExpression(cluster, verbose=FALSE)$RNA)
avg.exp$gene <- rownames(avg.exp)
avg.exp.male <- filter(avg.exp, gene %in% filtered.MaleUP)
avg.exp.female <- filter(avg.exp, gene %in% filtered.FemaleUP)
avg.exp.unchanged <- filter(avg.exp, gene %in% unchanged2)
pdf(file=paste0(output,i,"_scatterplot.pdf"))
p1 <- ggplot(avg.exp.male, aes(Male, Female)) + geom_point(color="blue") + geom_point(data=avg.exp.female, color="red") + geom_point(data=avg.exp.unchanged, alpha="0.01")
print(p1)
graphics.off()
})
}
