#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_PrimedvUnprimed_doubletsremoved/DiffPlots/UpPrimed/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_PrimedvUnprimed_doubletsremoved/PrimedvUnprimed_doubletsremoved_filtered1.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("Genelists/POA_UpPrimedvUnprimed.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.state <- paste(Idents(mySeurat), mySeurat$state, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.state"
for (g in genelist){
try(
dir.create(paste0(output,g)) 
)
for (i in 0:35){
try({
out.dir <- (paste0(output,"/",g,"/"))
pdf(paste0(out.dir,i,".pdf"))
ident1 <- paste0(i,"_Primed")
ident2 <- paste0(i,"_Unprimed")
idents <- c(ident1,ident2)
vln <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=1, combine=TRUE)
print(vln)
graphics.off()
}
)
}}

