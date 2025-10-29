#!/usr/bin/env Rscript

sampleID <- "MaleBNST"
directory <- "MaleBNSTIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_BusMergeTest/BNST_BusMergeTest_MvP_Hailmary_cDNAonly"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedBNST/cDNA_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=200)
Intact.data <- read_count_output("Kallisto_Bus_Output/IntactBNST/cDNA_counts_em/", name="output", tcc=FALSE)
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactMale", min.cells=3, min.features=200)
Intact$sex <- "Male"
Primed$sex <- "Female"
Intact$status <- "IntactMale"
Primed$status <- "PrimedFemale"
mySeurat <- merge(Intact, y=c(Primed),add.cell.ids=c("IntactMale","PrimedFemale"), project="Merged")
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress="orig.ident")
mySeurat <- RunPCA(mySeurat, feature=VariableFeatures(object=mySeurat))
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat,reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_statussplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="status")
dev.off()
pdf(paste0(output,"_UMAP_statuslabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="status")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

