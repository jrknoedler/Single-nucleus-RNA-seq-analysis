#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_BusMergeTest/POA_BusMergeTest"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedPOA/combined_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=1000)
Intact.data <- read_count_output("Kallisto_Bus_Output/IntactPOA/combined_counts_em/", name="output", tcc=FALSE)
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactMale", min.cells=3, min.features=1000)
Unprimed.data <- read_count_output("Kallisto_Bus_Output/UnprimedPOA/combined_counts_em/", name="output", tcc=FALSE)
Unprimed <- CreateSeuratObject(counts = Unprimed.data, project = "UnprimedFemale", min.cells=3, min.features=1000)
Castrate.data <- read_count_output("Kallisto_Bus_Output/CastratePOA/combined_counts_em/", name="output", tcc=FALSE)
Castrate <- CreateSeuratObject(counts = Castrate.data, project = "CastrateMale", min.cells=3, min.features=1000)
Intact$sex <- "Male"
Primed$sex <- "Female"
Castrate$sex <- "Male"
Unprimed$sex <- "Female"
Intact$status <- "IntactMale"
Primed$status <- "PrimedFemale"
Castrate$status <- "CastrateMale"
Unprimed$status <- "UnprimedFemale"
Intact$hormone <- "Steroidal"
Primed$hormone <- "Steroidal"
Castrate$hormone <- "Nonsteroidal"
Unprimed$hormone <- "Nonsteroidal"
mySeurat <- merge(Intact, y=c(Primed, Unprimed, Castrate),add.cell.ids=c("IntactMale","PrimedFemale","UnprimedFemale","CastrateMale"), project="Merged")
mySeurat <- subset(mySeurat, subset=nCount_RNA < 60000)
mySeurat
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress="orig.ident")
mySeurat <- RunPCA(mySeurat, feature=VariableFeatures(object=mySeurat))
warnings()
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.5)
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
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="status")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="hormone")
dev.off()
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))

