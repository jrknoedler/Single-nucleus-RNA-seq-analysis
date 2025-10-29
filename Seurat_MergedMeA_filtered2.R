#!/usr/bin/env Rscript

output <- "Seurat/MeA_IndependentAnalysis/MeA_MergedFinal2_Filtered2"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
#library(DoubletFinder)

mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_MergedFinal2_Filtered1.rds")
mySeurat <- subset(mySeurat, idents=c("9","17"), invert=TRUE)
mySeurat <- SCTransform(mySeurat, vars.to.regress="orig.ident")
mySeurat <- RunPCA(mySeurat, feature=VariableFeatures(object=mySeurat))
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.5)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
#sweep.res.list <- paramSweep_v3(mySeurat, PCs=1:30, sct=TRUE)
#sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
#bcmvn <- find.pK(sweep.stats)
#bcmvn
#nExp_poi <- round(0.15*nrow(mySeurat@meta.data))
#mySeurat <- doubletFinder_v3(mySeurat, PCs=1:30, pN=0.25, pK=0.14, nExp=nExp_poi, reuse.pANN=FALSE, sct=TRUE)
#head(mySeurat[[]])
saveRDS(mySeurat, file=(paste0(output, ".rds")))
#pdf(file=paste0(output,"doublettest.pdf"))
#DimPlot(mySeurat, group.by="DF.classifications_0.25_0.14_1476")
#dev.off()


