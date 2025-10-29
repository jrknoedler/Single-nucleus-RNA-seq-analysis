#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/Paper1stpass/20pcs/POA_IntactvPrimed_Cluster_mast_allgenes_0.2thresh"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test_20pcs.rds")
pdf(file=paste0(output,"FeatureScatter.pdf"))
FeatureScatter(mySeurat, feature2="nFeature_RNA", feature1="nCount_RNA", group.by="Hormone")
dev.off()
pdf(file=paste0(output,"UMIVln.pdf"), width=25)
VlnPlot(mySeurat, features=c("nCount_RNA"), pt.size=0, split.by="Hormone")
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
head(mySeurat[[]])
genelist <- read.table("Genelists/POA_MvPAll.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
mySeurat$celltype.Hormone <- paste(Idents(mySeurat), mySeurat$Hormone, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.Hormone"
for (i in 0:32){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i,"_Primed")
sex.dimorphism <- FindMarkers(mySeurat, ident.1 = ident1, ident.2=ident2, min.pct=0.1, logfc.threshold=0.2,  verbose=TRUE, test.use="MAST")
write.csv(sex.dimorphism, file=paste0(output,i,"_dimorphism.csv"))
})
}
