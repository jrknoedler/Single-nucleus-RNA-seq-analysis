#!/usr/bin/env Rscript

sampleID <- "MaleBNST"
directory <- "MaleBNSTIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_Prepaperrecluster_2"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedBNST/combined_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=500)
Primed
Intact.data <- read_count_output("Kallisto_Bus_Output/IntactBNST/combined_counts_em/", name="output", tcc=FALSE)
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactMale", min.cells=3, min.features=500)
Intact
Unprimed.data <- read_count_output("Kallisto_Bus_Output/UnprimedBNST/combined_counts_em/", name="output", tcc=FALSE)
Unprimed <- CreateSeuratObject(counts = Unprimed.data, project = "UnprimedFemale", min.cells=3, min.features=500)
Unprimed

#cluster and save sample 1
Primed <- subset(Primed, subset=nCount_RNA < 60000)
Primed
Primed <- SCTransform(Primed, verbose=TRUE)
Primed[["percent.Malat1"]] <- PercentageFeatureSet(Primed, pattern = "Malat1")
Primed <- SCTransform(Primed, verbose=TRUE, vars.to.regress=c("percent.Malat1"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
hvg <- Primed@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
Primed <- RunPCA(Primed)
warnings()
Primed <- FindNeighbors(Primed, features=hvg.final)
Primed <- FindClusters(Primed, resolution=1.5)
Primed <- RunUMAP(Primed, dims=1:30)
pdf(paste0(output,"Gad1Primed.pdf"))
VlnPlot(Primed, features=c("Gad1"), pt.size=0)
dev.off()
pdf(paste0(output,"_Primed_UMAP.pdf"))
DimPlot(Primed, reduction="umap", label=TRUE)
dev.off()
Primed.markers <- FindAllMarkers(Primed, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Primed.markers, file=paste0(output,"_Primed_allposmarkers.csv"))
pdf(paste0(output,"_Primed_TopMarkerheatmap.pdf"))
top10 <- Primed.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Primed, features=top10$gene) + NoLegend()
dev.off()
pdf(paste0(output,"_primedUMIvln.pdf"), width=30)
VlnPlot(Primed,features=c("nCount_RNA"), pt.size=1)
dev.off() 
saveRDS(Primed, file=(paste0(output, "_Primed.rds")))

#cluster and save sample 2
Unprimed <- subset(Unprimed, subset=nCount_RNA < 60000)
Unprimed <- SCTransform(Unprimed, verbose=TRUE)

Unprimed[["percent.Malat1"]] <- PercentageFeatureSet(Unprimed, pattern = "Malat1")
Unprimed <- SCTransform(Unprimed, verbose=TRUE, vars.to.regress=c("percent.Malat1"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
hvg <- Unprimed@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
Unprimed <- RunPCA(Unprimed)

Unprimed <- RunPCA(Unprimed, features=hvg.final)
warnings()
Unprimed <- FindNeighbors(Unprimed, dims=1:30)
Unprimed <- FindClusters(Unprimed, resolution=1.5)
Unprimed <- RunUMAP(Unprimed, dims=1:30)
pdf(paste0(output,"_Unprimed_UMAP.pdf"))
DimPlot(Unprimed, reduction="umap", label=TRUE)
dev.off()
Unprimed.markers <- FindAllMarkers(Unprimed, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Unprimed.markers, file=paste0(output,"_Unprimed_allposmarkers.csv"))
pdf(paste0(output,"_Unprimed_TopMarkerheatmap.pdf"))
top10 <- Unprimed.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Unprimed, features=top10$gene) + NoLegend()
dev.off()
pdf(paste0(output,"_unprimedUMIvln.pdf"), width=30)
VlnPlot(Unprimed,features=c("nCount_RNA"), pt.size=1)
dev.off() 
saveRDS(Unprimed, file=(paste0(output, "_Unprimed.rds")))

#cluster and save sample 3
Intact <- subset(Intact, subset=nCount_RNA < 60000)
Intact <- SCTransform(Intact, verbose=TRUE)

Intact[["percent.Malat1"]] <- PercentageFeatureSet(Intact, pattern = "Malat1")
Intact <- SCTransform(Intact, verbose=TRUE, vars.to.regress=c("percent.Malat1"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
hvg <- Intact@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
Intact <- RunPCA(Intact, features=hvg.final)

warnings()
Intact <- FindNeighbors(Intact, dims=1:30)
Intact <- FindClusters(Intact, resolution=1.5)
Intact <- RunUMAP(Intact, dims=1:30)
pdf(paste0(output,"_Intact_UMAP.pdf"))
DimPlot(Intact, reduction="umap", label=TRUE)
dev.off()
Intact.markers <- FindAllMarkers(Intact, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(Intact.markers, file=paste0(output,"_Intact_allposmarkers.csv"))
pdf(paste0(output,"_Intact_TopMarkerheatmap.pdf"))
top10 <- Intact.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(Intact, features=top10$gene) + NoLegend()
dev.off()
pdf(paste0(output,"_maleUMIvln.pdf"), width=30)
VlnPlot(Intact,features=c("nCount_RNA"), pt.size=1)
dev.off() 
saveRDS(Intact, file=paste0(output, "_Intact.rds"))

