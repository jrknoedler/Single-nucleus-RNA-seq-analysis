#!/usr/bin/env Rscript


output <- "Seurat/POA_IndependentAnalysis/POA_Primed_Kiss1opt_500"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)

spliced <- read_count_output("Kallisto_Bus_Kiss1Opt_Out/POA_Primed/counts_unfiltered/", name="spliced", tcc=FALSE)
unspliced <- read_count_output("Kallisto_Bus_Kiss1Opt_Out/POA_Primed/counts_unfiltered/", name="unspliced", tcc=FALSE)

common_cells <- intersect(colnames(spliced), colnames(unspliced))
spliced <- spliced[, colnames(spliced) %in% common_cells]
unspliced <- unspliced[, colnames(unspliced) %in% common_cells]

all.equal(colnames(spliced), colnames(unspliced))
all.equal(rownames(spliced), rownames(unspliced))

Merged <- spliced + unspliced

rownames(Merged) <- rownames(spliced)
colnames(Merged) <- colnames(spliced)

Intact <- CreateSeuratObject(counts = Merged, project = "IntactMale", min.cells=3, min.features=500)
Intact
head(Intact[[]])
#pdf(file=paste0(output,"totalEsr1.pdf"))
#lnPlot(Intact, features=c("ENSMUSG00000019768.16"), pt.size=1)
#dev.off()
#pdf(file=paste0(output,"totalKiss1.pdf"))
#VlnPlot(Intact, features=c("ENSMUSG00000116158.1"),  pt.size=1)
#dev.off()
Intact <- SCTransform(Intact, verbose=TRUE)
Intact <- RunPCA(Intact)
warnings()
Intact <- FindNeighbors(Intact, dims=1:30)
Intact <- FindClusters(Intact, resolution=1)
Intact <- RunUMAP(Intact, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(Intact, reduction="umap", label=TRUE)
dev.off()
Intact <- BuildClusterTree(Intact, dims=1:30)
pdf(paste0(output,"_Intact_clustertree.pdf"))
PlotClusterTree(object=Intact)
dev.off()
Intact.markers <- FindAllMarkers(Intact, only.pos=TRUE, min.pct=0.25, min.diff.pct=0.2, logfc.threshold = 0.25, test.use="MAST")
write.csv(Intact.markers, file=paste0(output,"_Intact_allposmarkers.csv"))
pdf(paste0(output,"__TopMarkerheatmap.pdf"))
top10 <- Intact.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Intact, features=top10$gene) + NoLegend()
dev.off()
saveRDS(Intact, file=(paste0(output, ".rds")))