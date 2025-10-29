#!/usr/bin/env Rscript

#sampleID <- "PregPOA"
output <- "Seurat/MeA_Barcode_Transfer/BNST_Barcode_Lift_SCT"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(Matrix)
Primed.data <- read_count_output("Updated_Ref_Kb_Out/MeA_Primed/output_mtx_em", name="output", tcc=FALSE)
dim(Primed.data)
Unprimed.data <- read_count_output("Updated_Ref_Kb_Out/MeA_Unprimed/output_mtx_em", name="output", tcc=FALSE)
dim(Unprimed.data)
Male.data <- read_count_output("Updated_Ref_Kb_Out/Male_MeA/output_mtx_em", name="output", tcc=FALSE)
dim(Male.data)
Primed.barcodes <- read.table("Barcodes_For_Transfer/Primed_MeA_Barcodes.txt", header=FALSE)
Primed.barcodes <- unlist(Primed.barcodes)
head(Primed.barcodes)
Unprimed.barcodes <- read.table("Barcodes_For_Transfer/New_Unprimed_MeA_Barcodes.txt", header=FALSE)
Unprimed.barcodes <- unlist(Unprimed.barcodes)
head(Unprimed.barcodes)
Male.barcodes <- read.table("Barcodes_For_Transfer/Male_MeA_Barcodes.txt", header=FALSE)
Male.barcodes <- unlist(Male.barcodes)
head(Male.barcodes)
Primed.filtered <- Primed.data[,Primed.barcodes]
dim(Primed.filtered)
Male.filtered <- Male.data[,Male.barcodes]
dim(Male.filtered)
#Unprimed.filtered <- Unprimed.data[,Unprimed.barcodes]
#dim(Unprimed.filtered)
PrimedSeurat <- CreateSeuratObject(counts = Primed.filtered, project = "PrimedMeA", min.cells=3, min.features=500)
PrimedSeurat[["percent.Malat1"]] <- PercentageFeatureSet(PrimedSeurat, pattern = "Malat1")
PrimedSeurat <- SCTransform(PrimedSeurat, verbose=TRUE, vars.to.regress=c("percent.Malat1"))
UnprimedSeurat <- CreateSeuratObject(counts = Unprimed.data, project = "UnprimedMeA", min.cells=3, min.features=500)
UnprimedSeurat <- subset(UnprimedSeurat, cells=Unprimed.barcodes)
UnprimedSeurat

UnprimedSeurat[["percent.Malat1"]] <- PercentageFeatureSet(UnprimedSeurat, pattern = "Malat1")
UnprimedSeurat <- SCTransform(UnprimedSeurat, verbose=TRUE, vars.to.regress=c("percent.Malat1"))
MaleSeurat <- CreateSeuratObject(counts=Male.filtered, project = "MaleMeA", min.cells=3, min.features=500)
MaleSeurat[["percent.Malat1"]] <- PercentageFeatureSet(MaleSeurat, pattern = "Malat1")
MaleSeurat <- SCTransform(MaleSeurat, verbose=TRUE, vars.to.regress=c("percent.Malat1"))
PrimedSeurat$sex <- "Female"
PrimedSeurat$Hormone <- "Primed"
UnprimedSeurat$sex <- "Female"
UnprimedSeurat$Hormone <- "Unprimed"
MaleSeurat$sex <- "Male"
MaleSeurat$Hormone <- "Intact"
Pregnant <- readRDS("Seurat/PregMeA_Exploratory/PregMeA_Filtered4.rds")
Pregnant$sex <- "Female"
Pregnant$Hormone <- "Pregnant"
Pregnant[["percent.Malat1"]] <- PercentageFeatureSet(Pregnant, pattern = "Malat1")

list <- c(MaleSeurat,UnprimedSeurat,PrimedSeurat,Pregnant)
Int.features <- SelectIntegrationFeatures(object.list=list, nfeatures=3000)
Sample.list <- PrepSCTIntegration(list, anchor.features=Int.features, verbose=TRUE)
Int.anchors <- FindIntegrationAnchors(Sample.list, normalization.method="SCT", anchor.features=Int.features, verbose=TRUE)
mySeurat <- IntegrateData(anchorset=Int.anchors, normalization.method="SCT", verbose=FALSE)
mySeurat <- RunPCA(mySeurat)


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
