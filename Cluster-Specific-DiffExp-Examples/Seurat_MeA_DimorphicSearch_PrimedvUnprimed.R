#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/Paper_DraftAnalysis/MAST/PvU/"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DESeq2)


mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")


data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)



datas <- SampleUMI(data, max.umi=40000, upsample=FALSE, verbose=TRUE)
datar <- round(datas)

dsSeurat <- CreateSeuratObject(datar, project="IntegerMatrix", assay="RNA")

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone


dsSeurat <- AddMetaData(dsSeurat, metadata=idents, col.name="Ident")
dsSeurat <- AddMetaData(dsSeurat, metadata=hormone, col.name="Hormone")
Idents(dsSeurat) <- "Ident"
head(dsSeurat[[]])
DefaultAssay(dsSeurat) <- "RNA"
dsSeurat <- NormalizeData(dsSeurat)
head(dsSeurat[[]])
genelist <- read.table("Genelists/MeA_PvU_1.5cutoff.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
dim(x=dsSeurat)
genes.10x <- (x=rownames(x=dsSeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
dsSeurat$celltype.Hormone <- paste(Idents(dsSeurat), mySeurat$Hormone, sep="_")
dsSeurat$celltype <- Idents(dsSeurat)
Idents(dsSeurat) <- "celltype.Hormone"
total_dimorphic = data.frame()
for (i in 0:37){
try({
dsSeurat <- subset(dsSeurat, nCount_RNA > 7500)
ident1 <- paste0(i,"_Primed")
ident2 <- paste0(i,"_Unprimed")
sex.dimorphism <- FindMarkers(dsSeurat, assay="RNA", ident.1 = ident1, ident.2=ident2, features=filtered.genelist, logfc.threshold=0, min.pct=0, only.pos=FALSE, test.use="MAST", verbose=TRUE)
write.csv(sex.dimorphism, file=paste0(output,i,"_dimorphism.csv"))
sex.dimorphism <- data.frame(sex.dimorphism)
#eDEGs <- sex.dimorphism[rownames(sex.dimorphism) %in% filtered.genelist,]
#rawpcounts <- nrow(eDEGs[eDEGs$p_val < 0.05,])
padjcounts <- nrow(eDEGs[eDEGs$p_val_adj < 0.05,])
df <- data.frame(i, padjcounts)
total_dimorphic=rbind(total_dimorphic, df)
})
}
write.csv(total_dimorphic, file=paste0(output,"_TRAPsDEGsdibyclust.csv"))
