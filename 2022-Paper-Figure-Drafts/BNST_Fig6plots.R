#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/BNST_IndependentAnalysis/BNST_Plotsupdate_regress30pcres1.2"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered2_sexclude_malat1regress_30pcsres1.2FINAL.rds")
pdf(file=paste0(output,"Malat1counts_Vln.pdf"), width=40)
VlnPlot(mySeurat, features=c("Malat1"), pt.size=0, slot="data", split.by="Hormone")
dev.off()
DefaultAssay(mySeurat) <- "RNA"
genelist <- read.table("Genelists/BNST_MvP_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist

mySeurat <- NormalizeData(mySeurat)
Avg <- AverageExpression(mySeurat, assays=c("RNA"), features=filtered.genelist)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/MvPsDEGclusteravgs.csv")
pdf(file=paste0(output,"Cyp19a1vlntest_dimorphicclusters.pdf"), width=24)
VlnPlot(mySeurat, features=c("Cyp19a1"), pt.size=0, split.by="Hormone", idents=c("3","10","13","14","16","17"), cols=c("light blue","pink", "light green"))
dev.off()
Aro.markers <- read.table("/scratch/users/knoedler/DotPlotMarkerLists/BNST_AroMarkers.txt")
Aro.markers <- unlist(Aro.markers)
pdf(file=paste0(output,"Cyp19a1_MarkerVln.pdf"), height=18)
VlnPlot(mySeurat, features=Aro.markers, idents=c("3","10","13","14","16","17","25"), ncol=1, pt.size=0)
dev.off()
props <- prop.table(table(Idents(mySeurat), mySeurat$Hormone), margin=2)
write.csv(props, file=paste0(output, "BNST_propperclust.csv"))
counts <- table(Idents(mySeurat), mySeurat$Hormone)
write.csv(counts, file=paste0(output,"BNST_countsperclust.csv"))
pdf(file=paste0(output,"Clusterbias.pdf"))
DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="blue", "4"="blue", "5"="light gray", "6"="light gray", 
	"7"="light gray", "8"="light gray", "9"="light gray", "10"="light green", "11"="light gray", "12"="light gray", "13"="light gray", 
	"14"="light gray", "15"="light gray", "16"="light gray", "17"="blue", "18"="light gray", "19"="light gray", "20"="orange", 
	"21"="light gray", "22"="light gray", "23"="pink", "24"="light gray", "25"="light gray", "26"="light gray", "27"="light gray", 
	"28"="light gray", "29"="light gray", "30"="light gray", "31"="light gray", "32"="light gray", "33"="light gray","34"="light gray", "35"="light gray"), label=FALSE, reduction="umap")
dev.off() 
#pdf(file=paste0(output,"Primedcluster.pdf"))
#DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="light gray", "4"="light gray", "5"="light gray", "6"="light gray", 
#	"7"="light gray", "8"="light gray", "9"="light gray", "10"="light gray", "11"="light gray", "12"="light gray", "13"="Red", 
#	"14"="light gray", "15"="light gray", "16"="light gray", "17"="light gray", "18"="light gray", "19"="light gray", "20"="light gray", 
#	"21"="light gray", "22"="Red", "23"="light gray", "24"="light gray", "25"="light gray", "26"="light gray", "27"="light gray", 
#	"28"="light gray", "29"="light gray", "30"="light gray", "31"="light gray", "32"="light gray", "33"="light gray"), label=FALSE, reduction="umap")
#dev.off() 
#pdf(file=paste0(output,"Femaleclusters.pdf"))
#DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="light gray", "4"="orange", "5"="light gray", "6"="light gray", 
#	"7"="orange", "8"="light gray", "9"="light gray", "10"="light gray", "11"="light gray", "12"="light gray", "13"="light gray", 
#	"14"="light gray", "15"="light gray", "16"="light gray", "17"="light gray", "18"="light gray", "19"="light gray", "20"="light gray", 
#	"21"="light gray", "22"="light gray", "23"="light gray", "24"="light gray", "25"="light gray", "26"="light gray", "27"="light gray", 
#	"28"="light gray", "29"="light gray", "30"="light gray", "31"="light gray", "32"="light gray", "33"="light gray"), label=FALSE, reduction="umap")
#dev.off() 

genelist <- read.table("Genelists/BNST_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
MvP <- read.table("Genelists/BNST_MvP_1.5cutoff.txt", header=FALSE)
MvU <- read.table("Genelists/BNST_MvU_1.5cutoff.txt", header=FALSE)
sDEGs <- rbind(MvP, MvU)
sDEGs <- unlist(sDEGs)
sDEGs <- unique(sDEGs)
filtered.sDEGs <- intersect(sDEGs, genes.10x)
all <- rbind(genelist, MvP, MvU)
all <- unlist(all)
all <- unique(all)
filtered.all <- intersect(all, genes.10x)
singlesig <- read.table("/scratch/users/knoedler/Heatmaps/Seurat/MetaDEGs/BNST_SingleDimorphic.txt", header=FALSE)
singlesig <- unlist(singlesig)
singlesig.filtered <- setdiff(filtered.all, singlesig)
mySeurat <- NormalizeData(mySeurat)
sex.markers <- FindAllMarkers(mySeurat, features=filtered.all, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST", return.thresh=1e-4)
write.csv(sex.markers, file=paste0(output,"Sexmarkerstats.csv"))
mySeurat2 <- BuildClusterTree(mySeurat, reorder=TRUE, dims=1:30)

ScaledSeurat <- ScaleData(mySeurat2, features=rownames(mySeurat), do.center=TRUE, do.scale=TRUE)

pdf(paste0(output,"_SexMarkerheatmapreordered.pdf"))
top10 <- sex.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(ScaledSeurat, features=top10$gene, raster=FALSE) 
dev.off()
top10genes <- top10$gene
top10genes <- unique(top10genes)
Avg <- AverageExpression(mySeurat, assays=c("RNA"), features=top10genes)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/BNST/AllSigDEGclusteravgs.csv")

Avg2 <- AverageExpression(mySeurat, assays=c("RNA"), features=singlesig.filtered)
write.csv(Avg2, file="/scratch/users/knoedler/Heatmaps/Seurat/MetaDEGs/BNSTMetaDEGs.csv")

pdf(file=paste0(output,"Malat1_Vln.pdf"), width=40)
VlnPlot(mySeurat, features=c("Malat1"), pt.size=0, split.by="Hormone")
dev.off()

pdf(file=paste0(output,"Tac1clust14_Vln.pdf"))
VlnPlot(mySeurat, features=c("Tac1"), idents=c("14"), pt.size=0, split.by="Hormone", cols=c("light blue","pink", "light green"))
dev.off()
pdf(file=paste0(output,"Tac1combined.pdf"), width=25)
VlnPlot(mySeurat, features=c("Tac1"), pt.size=0, split.by="Hormone")
dev.off()
pdf(file=paste0(output,"Tac1pooled.pdf"), width=25)
VlnPlot(mySeurat, features=c("Tac1"), pt.size=0)
dev.off()
pdf(file=paste0(output,"Tac1raw.pdf"), width=25)
VlnPlot(mySeurat, features=c("Tac1"), pt.size=1)
dev.off()
pdf(file=paste0(output,"Cyp19a1_Vln.pdf"), width=40)
VlnPlot(mySeurat, features=c("Cyp19a1"), pt.size=0, split.by="Hormone")
dev.off()
pdf(file=paste0(output, "Slc17a6_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Slc17a6"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Cyp19a1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Cyp19a1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Tac1_featureplot2.pdf"))
FeaturePlot(mySeurat, features=c("Tac1"), order=FALSE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Hormonelabel_UMAP_featureplot.pdf"))
DimPlot(mySeurat, group.by="Hormone", label=FALSE, reduction="umap")
dev.off()
pdf(file=paste0(output, "UMAP_nolabel.pdf"))
DimPlot(mySeurat, label=FALSE, reduction="umap")
dev.off()
pdf(file=paste0(output, "Slc32a1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Slc32a1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Slc18a2_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Slc18a2"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Nrip1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Nrip1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Gad1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Gad1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Pappa_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Pappa"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "St18_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("St18"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Col25a1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Col25a1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Esr2_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Esr2"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Esr2_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Esr2"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Th_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Th"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Nptx2_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Nptx2"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Soc2Vln_featureplot.pdf"))
VlnPlot(mySeurat, features=c("Socs2"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
dev.off()
Socs2 <- VlnPlot(mySeurat, features=c("Socs2"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
sink(file=paste0(output,"Socs2counts.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Socs2$data
sink()

pdf(file=paste0(output, "TroVln_featureplot.pdf"))
VlnPlot(mySeurat, features=c("Tro"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
dev.off()
Tro <- VlnPlot(mySeurat, features=c("Tro"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
sink(file=paste0(output,"Trocounts.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Tro$data
sink()

pdf(file=paste0(output, "Klhl1Vln_featureplot.pdf"))
VlnPlot(mySeurat, features=c("Klhl1"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
dev.off()
Klhl1 <- VlnPlot(mySeurat, features=c("Klhl1"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
sink(file=paste0(output,"Klhl1counts.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Klhl1$data
sink()

pdf(file=paste0(output, "Klhl1Vln_featureplot.pdf"))
VlnPlot(mySeurat, features=c("Klhl1"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
dev.off()
Klhl1 <- VlnPlot(mySeurat, features=c("Klhl1"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
sink(file=paste0(output,"Klhl1counts.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Klhl1$data
sink()

pdf(file=paste0(output, "Mical2Vln_featureplot.pdf"))
VlnPlot(mySeurat, features=c("Mical2"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
dev.off()
Mical2 <- VlnPlot(mySeurat, features=c("Mical2"), order=TRUE, idents=c("2"), split.by="Hormone", pt.size=0)
sink(file=paste0(output,"Mical2counts.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Mical2$data
sink()


excitatory.use <- WhichCells(mySeurat, idents=c("19","30","7","22"))
mySeurat <- SetIdent(mySeurat, cells=excitatory.use, value="Excitatory")
idents <- mySeurat@active.ident
idents
inhibitory.use <- WhichCells(mySeurat, idents=c("1","16","31","3","2","18","23","5","14","8","17","25","27","13","4","15","9","6","0","26","29","12","10","21","20","11","28","24","32"))
mySeurat <- SetIdent(mySeurat, cells=inhibitory.use, value="Inhibitory")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
seDEGs <- read.table("RNASeqKmeans/RNASeqKmeans/BNST_1.5unique.txt")
unlist(seDEGs)
seDEGs <- as.matrix(seDEGs)
seDEGs.filtered <- intersect(seDEGs, genes.10x)
DiffMarks <- FindMarkers(mySeurat, ident.1="Excitatory", ident.2 = "Inhibitory", min.pct=0, only.pos=FALSE, logfc.threshold=0, test.use="MAST", features=seDEGs.filtered)
write.csv(DiffMarks, file=paste0(output,"_ExcitatoryvInhibitory.csv"))