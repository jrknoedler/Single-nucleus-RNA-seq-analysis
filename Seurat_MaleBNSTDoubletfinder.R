output <- "Seurat/MaleBNST_Preprocessing/MaleBNST_Doubletfinder_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
mySeurat <- readRDS("Seurat/MaleBNST_Preprocessing/MaleBNST_filtered_redo.rds")
mySeurat
sweep.res.list_Male <- paramSweep_v3(mySeurat, PCs=1:30, sct=TRUE)
sweep.stats_Male <- summarizeSweep(sweep.res.list_Male, GT=FALSE)
bcmvn_Male <- find.pK(sweep.stats_Male)
annotations <- mySeurat@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.01*2700)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
seu_Male <- doubletFinder_v3(mySeurat, PCs=1:30, pN=0.25, pK= 0.09, nExp=nExp_poi , reuse.pANN=FALSE, sct=TRUE)
seu_Male <- doubletFinder_v3(mySeurat, PCs=1:30, pN=0.25, pK= 0.09, nExp=nExp_poi.adj , reuse.pANN=FALSE, sct=TRUE)
head(seu_Male[[]])
saveRDS(seu_Male, paste0(output,"doublets.RDS"))
