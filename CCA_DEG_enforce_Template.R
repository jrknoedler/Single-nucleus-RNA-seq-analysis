##This script merges two Seurat object with canonical correlation analysis, with the stipulation that any genes that are DEGs between conditions of interest
##are manually added to the list of integration features. CCA automatically excludes genes that are present in one dataset but not the other, which
##can overfit datasets and inhibit identification of condition-specific cell types. The exact lists of DEGs will depend on which datasets you're 
##integrating. 


##Specific output folder
output <- "Path/to/output/prefix_"

##Load libraries
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

##Load dataset 1
Published <- readRDS("Path/to/Published.rds")

##Load dataset 2
Pregnant <- readRDS("Path/to/Pregnant.rds")

##Create list of objects to be integrated
list <- c(Published, Pregnant)

##Get all genes present in Published data
genes.published <- (x=rownames(x=Published))

##Get all genes present in Pregnant data
genes.pregnant <- (x=rownames(x=Pregnant))

##Get only genes present in BOTH datasets
allgenes <- intersect(genes.published, genes.pregnant)

##Get relevant list of DEGs from TRAPseq experiment (e.g., all genes that are DEGs vs pregnant animals). You can also combine multiple lists here.
##Note that you may need to change the command based on file format (e.g. "read.table" vs "read.cvs"), and check headers etc to make
##sure you get everything.
DEGs <- read.csv("path/to/genes.txt", header=FALSE)

##reformat object
DEGs <- unlist(DEGs)

##Remove any duplicates on list (like if you combined four partially overlapping lists)
DEGs <- unique(DEGs)

##Get all DEGs that are detected in snRNAseq datset. 
DEGs.filtered <- intersect(DEGs, allgenes)

##Get genes to exclude from clustering (Xist, Tsix and y-chromosome genes)
genelist <- read.table("path/to/ExcludeIDs.txt", header=FALSE)
genelist <- unlist(genelist)
genelist <- as.matrix(genelist)



##This particular script uses standard z-scaling rather than SCT.
##Basically, Abbas Rizvi was sure somebody would ask us to integrate by CCA
##rather than naive integration, or at least show that we could still find the
##same or similar cell types (especially so-called 'estrus specific' types
##like VMHvl Esr1/Cckar). It turned out that z-scaling worked better for this purpose.
##Feel free to try it this way or use the SCT integration workflow.

#Set default assay to RNA instead of SCT
DefaultAssay(Published) <- "RNA"
DefaultAssay(Pregnant) <- "RNA"

#Normalize data
Published <- NormalizeData(Published)
Pregnant <- NormalizeData(Pregnant)

Published <- FindVariableFeatures(Published, selection.method="vst", nfeatures=2000)

UnprimedPOA <- FindVariableFeatures(UnprimedPOA, selection.method="vst", nfeatures=2000)

MalePOA <- FindVariableFeatures(MalePOA, selection.method="vst", nfeatures=2000)

features <- SelectIntegrationFeatures(object.list=list)

features.final <- union(features, DEGs.filtered)
features.final < unlist(features.final)
POA.anchors <- FindIntegrationAnchors(object.list=list, anchor.features = features)

POA.comb <- IntegrateData(anchorset=POA.anchors)

DefaultAssay(POA.comb) <- "integrated"

POA.comb <- ScaleData(POA.comb)
POA.comb <- RunPCA(POA.comb)
POA.comb <- RunUMAP(POA.comb, reduction="pca", dims=1:30)
POA.comb <- FindNeighbors(POA.comb, reduction="pca", dims=1:30)
POA.comb <- FindClusters(POA.comb, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(POA.comb, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
DefaultAssay(POA.comb) <- "RNA"
POA.comb <- NormalizeData(POA.comb)
pdf(file=paste0(output,"MarkerVln.pdf"), width=32, height=20)
VlnPlot(POA.comb, features=c("Esr1","Pgr","Slc17a6","Gad1","Kiss1"), pt.size=0, ncol=1)
dev.off()
pdf(paste0(output,"_Kiss1vln.pdf"), width=20)
VlnPlot(POA.comb, features=c("Kiss1"), split.by="sex", pt.size=0)
dev.off()
mySeurat.markers <- FindAllMarkers(POA.comb, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(POA.comb, file=paste0(output,"CCAscaledata_DEGsonly.rds"))