#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/Vln_Feb/"
MvPOutput <- "Seurat/POA_IndependentAnalysis/Vln_Feb/MvP/"
PvUOutput <- "Seurat/POA_IndependentAnalysis/Vln_Feb/PvU/"
MvUOutput <- "Seurat/POA_IndependentAnalysis/Vln_Feb/MvU/"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat
mySeurat <- NormalizeData(mySeurat)
genelist <- read.table("topGO/Total_ByRegion/POA_Genesonly.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
mySeurat$celltype.Hormone <- paste(Idents(mySeurat), mySeurat$Hormone, sep="_")
mySeurat$celltype <- Idents(mySeurat)
Idents(mySeurat) <- "celltype.Hormone"
for (g in genelist){
pdf(file=paste0(MvPOutput,g,".pdf"))
for (i in 0:38){
try({
ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i, "_Primed")
#ident3 <- paste0(i, "_Unprimed")
idents <- c(ident1,ident2)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=0, combine=TRUE, cols=c("light blue", "pink"))
plot(vlnvar)
})
}
graphics.off()
}

for (g in genelist){
pdf(file=paste0(MvUOutput,g,".pdf"))
for (i in 0:38){
try({
ident1 <- paste0(i,"_Intact")
#ident2 <- paste0(i, "_Primed")
ident3 <- paste0(i, "_Unprimed")
idents <- c(ident1,ident3)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=0, combine=TRUE, cols=c("light blue", "light green"))
plot(vlnvar)
})
}
graphics.off()
}

for (g in genelist){
pdf(file=paste0(PvUOutput,g,".pdf"))
for (i in 0:38){
try({
#ident1 <- paste0(i,"_Intact")
ident2 <- paste0(i, "_Primed")
ident3 <- paste0(i, "_Unprimed")
idents <- c(ident2,ident3)
vlnvar <- paste0(i,"vln")
vlnvar <- VlnPlot(mySeurat, features=c(g), idents=idents , pt.size=0, combine=TRUE, cols=c("pink", "light green"))
plot(vlnvar)
})
}
graphics.off()
}
