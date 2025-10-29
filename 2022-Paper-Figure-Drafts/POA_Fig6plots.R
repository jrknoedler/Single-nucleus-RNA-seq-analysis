#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_Kiss1Opt"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(file=paste0(output,"Malecluster.pdf"))
DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="light gray", "4"="light gray", "5"="light gray", "6"="light gray", 
	"7"="light gray", "8"="pink", "9"="light gray", "10"="blue", "11"="blue", "12"="orange", "13"="blue", 
	"14"="light gray", "15"="light gray", "16"="light gray", "17"="light gray", "18"="light gray", "19"="blue", "20"="blue", 
	"21"="blue", "22"="light gray", "23"="light gray", "24"="blue", "25"="orange", "26"="orange", "27"="light gray", 
	"28"="light gray", "29"="light gray", "30"="light gray", "31"="light gray", "32"="light gray", "33"="blue", "34"="light gray", "35"="light gray","36"="light gray", "37"="light gray", "38"="light gray"), label=FALSE, reduction="umap")
dev.off() 
#pdf(file=paste0(output,"Primedcluster.pdf"))
#DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="light gray", "4"="light gray", "5"="light gray", "6"="light gray", 
	"7"="light gray", "8"="light gray", "9"="Red", "10"="light gray", "11"="light gray", "12"="light gray", "13"="light gray", 
	"14"="light gray", "15"="light gray", "16"="light gray", "17"="light gray", "18"="light gray", "19"="light gray", "20"="light gray", 
	"21"="light gray", "22"="light gray", "23"="light gray", "24"="light gray", "25"="light gray", "26"="light gray", "27"="light gray", 
	"28"="light gray", "29"="light gray", "30"="light gray", "31"="light gray", "32"="light gray", "33"="light gray", "34"="light gray"), label=FALSE, reduction="umap")
#dev.off() 
#pdf(file=paste0(output,"Femaleclusters.pdf"))
#DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="light gray", "4"="orange", "5"="light gray", "6"="light gray", 
	"7"="light gray", "8"="light gray", "9"="light gray", "10"="light gray", "11"="light gray", "12"="light gray", "13"="light gray", 
	"14"="light gray", "15"="light gray", "16"="light gray", "17"="light gray", "18"="light gray", "19"="light gray", "20"="Orange", 
	"21"="Orange", "22"="light gray", "23"="light gray", "24"="light gray", "25"="light gray", "26"="light gray", "27"="light gray", 
	"28"="light gray", "29"="light gray", "30"="light gray", "31"="light gray", "32"="light gray", "33"="light gray", "34"="light gray"), label=FALSE, reduction="umap")
#dev.off() 
#pdf(file=paste0(output,"Clusterhighlights.pdf"))
#DimPlot(mySeurat, cols=c("0"="light gray", "1"="light gray", "2"="light gray", "3"="light gray", "4"="orange", "5"="light gray", "6"="light gray", 
	"7"="light gray", "8"="light gray", "9"="Red", "10"="light gray", "11"="light gray", "12"="light gray", "13"="light gray", 
	"14"="light gray", "15"="light gray", "16"="dark blue", "17"="light gray", "18"="light gray", "19"="light gray", "20"="Orange", 
	"21"="Orange", "22"="light gray", "23"="light gray", "24"="light gray", "25"="light gray", "26"="light gray", "27"="light gray", 
	"28"="light gray", "29"="light gray", "30"="light gray", "31"="light gray", "32"="light gray", "33"="light gray", "34"="light gray"), label=FALSE, reduction="umap")
#dev.off()
pdf(file=paste0(output, "Hormonelabel_UMAP_featureplot.pdf"))
DimPlot(mySeurat, group.by="Hormone", label=FALSE, reduction="umap")
dev.off()
pdf(file=paste0(output, "Slc32a1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Slc32a1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Slc18a2_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Slc18a2"), order=TRUE, cols=c("light blue", "red"))
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
excitatory.use <- WhichCells(mySeurat, idents=c("9","0","3","23","1","18","24","25","29","11"))
mySeurat <- SetIdent(mySeurat, cells=excitatory.use, value="Excitatory")
idents <- mySeurat@active.ident
idents
inhibitory.use <- WhichCells(mySeurat, idents=c("19","7","2","20","21","16","14","28","4","6","13","8","5","30","15","27","33"))
mySeurat <- SetIdent(mySeurat, cells=inhibitory.use, value="Inhibitory")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
seDEGs <- read.table("RNASeqKmeans/RNASeqKmeans/POA_1.5unique.txt")
unlist(seDEGs)
seDEGs <- as.matrix(seDEGs)
seDEGs.filtered <- intersect(seDEGs, genes.10x)
DiffMarks <- FindMarkers(mySeurat, ident.1="Excitatory", ident.2 = "Inhibitory", min.pct=0, only.pos=FALSE, logfc.threshold=0, test.use="MAST", features=seDEGs.filtered)
write.csv(DiffMarks, file=paste0(output,"_ExcitatoryvInhibitory.csv"))