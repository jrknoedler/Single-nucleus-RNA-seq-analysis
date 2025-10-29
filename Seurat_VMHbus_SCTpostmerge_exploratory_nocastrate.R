#!/usr/bin/env Rscript

sampleID <- "MaleVMH"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_BusMergeTest/VMH_BusMergeTest_nocastrate_fewPCs"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedVMH/combined_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=500)
Primed
Intact.data <- read_count_output("Kallisto_Bus_Output/IntactVMH/combined_counts_em/", name="output", tcc=FALSE)
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactMale", min.cells=3, min.features=500)
Intact
Unprimed.data <- read_count_output("Kallisto_Bus_Output/UnprimedVMH/combined_counts_em/", name="output", tcc=FALSE)
Unprimed <- CreateSeuratObject(counts = Unprimed.data, project = "UnprimedFemale", min.cells=3, min.features=500)
Unprimed
Intact$sex <- "Male"
Primed$sex <- "Female"
Unprimed$sex <- "Female"
Intact$status <- "IntactMale"
Primed$status <- "PrimedFemale"
Unprimed$status <- "UnprimedFemale"
Intact$hormone <- "Steroidal"
Primed$hormone <- "Steroidal"
Unprimed$hormone <- "Nonsteroidal"
mySeurat <- merge(Intact, y=c(Primed, Unprimed),add.cell.ids=c("IntactMale","PrimedFemale","UnprimedFemale"), project="Merged")
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress="orig.ident")
genelist <- read.table("Genelists/Sexchrom.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
mySeurat <- RunPCA(mySeurat, features=hvg.final)
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:10)
mySeurat <- FindClusters(mySeurat, resolution=1)
mySeurat <- RunUMAP(mySeurat, dims=1:10)
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
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="status")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="hormone")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

