output <- "Seurat/MalePOA_Preprocessing/MalePOA_Doubletfinder2"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
mySeurat <- readRDS("Seurat/MaleBNST_Preprocessing/MaleBNST_Doubletfinder_doublets.RDS")
mySeurat
UpdateSeuratObject(mySeurat)
head(mySeurat[[]])
pdf(file=paste0(output,"_doublets.pdf"))
DimPlot(mySeurat, reduction="umap", group.by="DF.classifications_0.25_0.09_27")
dev.off()

