#!/usr/bin/env Rscript

#sampleID <- "PregPOA"
output <- "Seurat/VMH_Barcode_Transfer/VMH_Barcode_Lift"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
Primed.data <- read_count_output("Updated_Ref_Kb_Out/VMH_Primed/output_mtx_em", name="output", tcc=FALSE)
Unprimed.data <- read_count_output("Updated_Ref_Kb_Out/VMH_Unprimed/output_mtx_em", name="output", tcc=FALSE)
Male.data <- read_count_output("Updated_Ref_Kb_Out/Male_VMH/output_mtx_em", name="output", tcc=FALSE)
Primed.barcodes <- read.table("Barcodes_For_Transfer/Primed_VMH_Barcodes.txt", header=FALSE)
Primed.barcodes <- unlist(Primed.barcodes)
Unprimed.barcodes <- read.table("Barcodes_For_Transfer/Unprimed_VMH_Barcodes.txt", header=FALSE)
Unprimed.barcodes <- unlist(Unprimed.barcodes)
Male.barcodes <- read.table("Barcodes_For_Transfer/Male_VMH_Barcodes.txt", header=FALSE)
Male.barcodes <- unlist(Male.barcodes)
Primed.filtered <- Primed.data[,Primed.barcodes]
Unprimed.filtered <- Unprimed.data[,Unprimed.barcodes]
Male.filtered <- Male.data[,Male.barcodes]
PrimedSeurat <- CreateSeuratObject(counts = Primed.filtered, project = "PrimedVMH", min.cells=3, min.features=500)
PrimedSeurat[["percent.Malat1"]] <- PercentageFeatureSet(PrimedSeurat, pattern = "Malat1")
head(PrimedSeurat[[]])
tail(PrimedSeurat[[]])
UnprimedSeurat <- CreateSeuratObject(counts = Unprimed.filtered, project = "UnprimedVMH", min.cells=3, min.features=500)
UnprimedSeurat[["percent.Malat1"]] <- PercentageFeatureSet(UnprimedSeurat, pattern = "Malat1")
head(UnprimedSeurat[[]])
tail(UnprimedSeurat[[]])
MaleSeurat <- CreateSeuratObject(counts=Male.filtered, project = "MaleVMH", min.cells=3, min.features=500)
MaleSeurat[["percent.Malat1"]] <- PercentageFeatureSet(MaleSeurat, pattern = "Malat1")
head(MaleSeurat[[]])
tail(MaleSeurat[[]])
PrimedSeurat$sex <- "Female"
PrimedSeurat$Hormone <- "Primed"
UnprimedSeurat$sex <- "Female"
UnprimedSeurat$Hormone <- "Unprimed"
MaleSeurat$sex <- "Male"
MaleSeurat$Hormone <- "Intact"
Pregnant <- readRDS("Seurat/PregVMH_Exploratory/PregVMH_Filtered7.rds")
Pregnant[["percent.Malat1"]] <- PercentageFeatureSet(Pregnant, pattern = "Malat1")
Pregnant$sex <- "Female"
Pregnant$Hormone <- "Pregnant"
mySeurat <- merge(x=Pregnant, y=c(PrimedSeurat,UnprimedSeurat,MaleSeurat), add.cell.ids=c("Pregnant","Primed","Unprimed","Male"), project="VMH_Preg_Merged")
#mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
head(mySeurat[[]])
tail(mySeurat[[]])
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.Malat1"))
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
mySeurat <- RunPCA(mySeurat, features=hvg.final)
mySeurat <- FindNeighbors(mySeurat, dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.2)
mySeurat <- RunUMAP(mySeurat, dims=1:30)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", group.by="Hormone", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_Sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", group.by="sex", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_HormoneSplit.pdf"), height=3.5, width=14)
DimPlot(mySeurat, reduction="umap", split.by="Hormone", label=TRUE)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
pdf(paste0(output,"RNAfeat.pdf"), height=7, width=42)
VlnPlot(mySeurat, features=("nCount_RNA"), pt.size=0)
dev.off()
pdf(paste0(output,"genenum.pdf"), height=7, width=42)
VlnPlot(mySeurat, features=c("nFeature_RNA"), pt.size=0)
dev.off()
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.30, logfc.threshold = 0.30, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
