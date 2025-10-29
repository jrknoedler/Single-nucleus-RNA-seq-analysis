#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/Paper_DraftAnalysis/MvP30pcs/BNST_MvP_1.5cutoffSCT_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")
#mySeurat <- subset(mySeurat, nCount_RNA >= 5000)

data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)



#datas <- SampleUMI(data, max.umi=15000, upsample=FALSE, verbose=TRUE)

#dsSeurat <- CreateSeuratObject(datas, project="IntegerMatrix", assay="RNA")

#idents <- mySeurat@active.ident
#hormone <- mySeurat$Hormone

#dsSeurat <- AddMetaData(dsSeurat, metadata=idents, col.name="Ident")
#dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
#Idents(mySeurat) <- "Ident"
#head(mySeurat[[]])
#DefaultAssay(mySeurat) <- "RNA"
#mySeurat <- NormalizeData(mySeurat)
#head(mySeurat[[]])
genelist <- read.table("Genelists/BNST_MvP_1.5cutoff.txt", header=FALSE)
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
for (i in 0:33){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i,"_Primed")
sex.dimorphism <- FindMarkers(mySeurat, assay="SCT", ident.1 = ident1, ident.2=ident2, features=filtered.genelist, min.cells.group=3, logfc.threshold=0.20, min.pct=0, only.pos=FALSE, verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(output,i,"_dimorphism.csv"))
})
}
