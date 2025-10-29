#!/usr/bin/env Rscript

sampleID <- "MaleVMH"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/UnprimedMeA_Exploratory/UnprimedMeA_emptyDrops_.01FDR"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Seurat.data <- read_count_output("Kallisto_Bus_Output/UnprimedMeA_Runscombined/combined_counts_em/", name="output", tcc=FALSE)
#!/usr/bin/env Rscript

sampleID <- "MaleVMH"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/UnprimedMeA_Exploratory/UnprimedMeA_1000cutoff_runscombined"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(DropletUtils)
mtx.data <- read_count_output("Kallisto_Bus_Output/UnprimedMeA_Runscombined/combined_counts_em/", name="output", tcc=FALSE)
set.seed(100)
e.out <- emptyDrops(mtx.data, lower=100)
e.out
is.cell <- e.out$FDR <=0.01
sum(is.cell, na.rm=TRUE)
table(Limited=e.out$Limited, Significant=is.cell)
pdf(file=paste0(output,"EmptyDropstest.pdf"))

plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")
dev.off()
mySeurat <- CreateSeuratObject(counts = is.cell, project = sampleID, min.cells=3, min.features=200)
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
mySeurat <- FindNeighbors(mySeurat, dims=1:20)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:20)
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
