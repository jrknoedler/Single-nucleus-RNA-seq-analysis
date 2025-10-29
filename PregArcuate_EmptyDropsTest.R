#!/usr/bin/env Rscript

sampleID <- "PregnantBNST"

output <- "Seurat/PregArc_Exploratory/PregArc_emptyDrops_750lower.001FDR_50000iter"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DropletUtils)

mtx.data <- read_count_output("Preg_KB_Out/Preg_Arc2/output_mtx_em/", name="output", tcc=FALSE)
head(mtx.data)
set.seed(100)
e.out <- emptyDrops(mtx.data, lower=750, niter=50000)
e.out
is.cell <- e.out$FDR <= 0.001
dim(is.cell)
sum(is.cell, na.rm=TRUE)
table(Limited=e.out$Limited, Significant=is.cell)
#saveRDS(is.cell, file="iscellint.rds")
#saveRDS(e.out, file="eoutint.rds")




filtered.counts <- mtx.data[,which(is.cell)]
dim(filtered.counts)
mySeurat <- CreateSeuratObject(counts = filtered.counts, project = sampleID, min.cells=3, min.features=500)
mySeurat
pdf(paste0(output,"_basic_counts.pdf"))
VlnPlot(mySeurat, features=c("nFeature_RNA", "nCount_RNA"))
dev.off()
pdf(paste0(output,"_Featurescatter_complexity.pdf"))
FeatureScatter(mySeurat, feature1="nCount_RNA", feature2="nFeature_RNA")
dev.off()
pdf(paste0(output,"_Esr1counts.pdf"))
VlnPlot(mySeurat, features = c("Esr1"))
dev.off()
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat <- SCTransform(mySeurat, verbose=TRUE)
mySeurat <- RunPCA(mySeurat, feature=VariableFeatures(object=mySeurat))
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output, "_umi.pdf"))
VlnPlot(mySeurat, features=c("nCount_RNA"), ncol=1, pt.size=0)
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))
