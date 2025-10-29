#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_MalevPrimed_SCTDE/Scatterplots/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(tibble)
mySeurat <- readRDS("Seurat/VMH_MalevPrimedMerged_Exploratory/VMH_MalevPrimed_Exploratory_filtered2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
head(mySeurat[[]])
mySeurat <- NormalizeData(mySeurat)
MaleUP <- read.table("Genelists/VMH_MalevPrimed_Maleup.txt", header=FALSE)
MaleUP <- unlist(MaleUP)
FemaleUP <- read.table("Genelists/VMH_MalevPrimed_Primedup.txt", header=FALSE)
FemaleUP <- unlist(FemaleUP)
for (i in 0:33){
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
avg.exp <- avg.exp %>% rownames_to_column('Gene')
avg.exp.male <- filter(avg.exp, gene %in% filtered.MaleUP)
avg.exp.male <- avg.exp.male %>% column_to_rownames('Gene')
avg.exp.male$category <- "Male"
avg.exp.female <- filter(avg.exp, gene %in% filtered.FemaleUP)
avg.exp.female <- avg.exp.female %>% column_to_rownames('Gene')
avg.exp.female$category <- "Female"
avg.exp.unchanged <- filter(avg.exp, gene %in% unchanged2)
avg.exp.unchanged <- avg.exp.unchanged %>% column_to_rownames('Gene') 
avg.exp.unchanged$category <- "Unchanged"
genes.to.label = c("Cckar","Phf21b","Setdb2","Pgr","Esr1","Mical2","Ezr","Pdyn","Rprm")
all.genes <- rbind(avg.exp.unchanged,avg.exp.male,avg.exp.female)
pdf(file=paste0(output,i,"_scatterplot.pdf"))
p1 <- ggplot(all.genes, aes(Male,Female)) + geom_point(aes(color=category, alpha=category))+scale_color_manual(values=c("Red","Blue","Black"))+scale_alpha_manual(values=c(1,1,0.05))
#p1 <- ggplot(avg.exp.female, aes(Male, Female)) + geom_point(color="red") + geom_point(data=avg.exp.male, color="blue") + geom_point(data=avg.exp.unchanged, alpha="0.01")
p1 <- LabelPoints(plot=p1, points=genes.to.label, repel=TRUE, xnudge=0, ynudge=0)
print(p1)
graphics.off()
})
}
