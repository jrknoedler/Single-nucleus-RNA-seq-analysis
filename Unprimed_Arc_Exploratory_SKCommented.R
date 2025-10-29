#!/usr/bin/env Rscript

### You'll need to make some modifications - this is a script to do a first pass on a single-cell RNAseq library. You'll need to look 
### at the output to make decisions about what to do with the data. Feel free to comment out any plots you think are unnecessary.

## set output directory
output <- "Seurat/Unprimed_Arcuate_Exploratory/Unprimed_Arcuate_1stPass_ParamsAdjusted"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)

### Read in the data. This is for BUStools output - you will need to use a different command (Read10X might work)
Seurat.data <- read_count_output("Preg_KB_Out/Unprimed_Arc/output_mtx_em/", name="output", tcc=FALSE)

### Turn count table into a Seurat Object
mySeurat <- CreateSeuratObject(counts = Seurat.data, project = sampleID, min.cells=3, min.features=500)

## Look at number of genes/cell before you cluster
pdf(paste0(output, "_basic_counts.pdf"))
VlnPlot(mySeurat, features=c("nFeature_RNA", "nCount_RNA"))
dev.off()

## Plot # of features vs # of genes per cell (can help identify doublets)
pdf(paste0(output, "_Featurescatter_complexity.pdf"))
FeatureScatter(mySeurat, feature1="nCount_RNA", feature2="nFeature_RNA")
dev.off()

## Look at what % of cells express Esr1 - you might want to replace this with Tacr1, or plot both)
pdf(paste0(output, "_Esr1counts.pdf"))
VlnPlot(mySeurat, features = c("Esr1"))
dev.off()

##Get rid of cells with too many counts (obvious doublets)
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)

##Run SCTransform (I use this but you can also just z-scale if you prefer)
mySeurat <- SCTransform(mySeurat, verbose=TRUE)

#######
##Sometimes you need to exclude genes, e.g. y-chromosome genes or other things you don't like, from PCA. I've commented out so you can ignore but you might ##need it someday 
##genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
##genelist <- unlist(genelist)
##genelist <- as.matrix(genelist)
##hvg <- mySeurat@assays$SCT@var.features
##hvg <- unlist(hvg)
##hvg <- as.matrix(hvg)
##hvg.final <- setdiff(hvg, genelist)
######

##Run PCA
mySeurat <- RunPCA(mySeurat)

##Run clustering. IMPORTANT: 30 PCs and resolution 1.5 is good when you're starting with thousands of cells but will almost
##certainly oversplit your PIPseq library. I'd recommend using fewer (maybe dims 1:10 and resolution 0.6) or testing a whole bunch

mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.5)
mySeurat <- RunUMAP(mySeurat, dims=1:30)

##make a UMAP
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()

###Plot some marker genes to see what types of cells you got
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=40)
VlnPlot(mySeurat, features=c("Gad1","Slc17a6","Esr1","Tacr1"), ncol=1, pt.size=0)
dev.off()

##Plot UMIs and genes by cluster
pdf(paste0(output,"RNAfeat.pdf"), height=7, width=42)
VlnPlot(mySeurat, features=("nCount_RNA","nFeature_RNA"), pt.size=0)
dev.off()

##Find marker genes with a generous threshold
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))

##Plot z-scored expression to visualize clusters that may be low-quality
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))
