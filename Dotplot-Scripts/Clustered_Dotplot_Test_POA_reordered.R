#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_CelltypesKiss1regressed_reordered"
regoutput <- "RegulonsCompiled/POA/Sexmarkersperclust/"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)


mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")
pdf(file=paste0(output,"plotnolab.pdf"))
DimPlot(mySeurat, label=FALSE, reduction="umap")
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
#pdf(file=paste0(output,"biasmarkers.pdf"), width=10, height=24)
#VlnPlot(mySeurat, idents=c("8","10","11","12","13","19","20","24","25","26","33"), features=rev(c("Gad1","Slc17a6","Angpt1","Ror1","Esr2","Ltbp2","Lmo7","Fbn2","Cdh23","Meis2","Klhl4","Kiss1","Rspo2")), ncol=1, pt.size=0)
#dev.off()
#pdf(file=paste0(output,"MaleClustMarkers.pdf"), width=10, height=30)
#VlnPlot(mySeurat, idents=c("10","19","24","33"), features=rev(c("Gad1","Slc17a6","Slc18a2","Cyp19a1","Shisal2b","Ror1","Fbn2","Meis2","Rspo2")), cols=c("dark green","dark blue","purple","brown"),ncol=1, pt.size=0)
#dev.off()
#pdf(file=paste0(output,"FemaleClustMarkers.pdf"), height=20)
#VlnPlot(mySeurat, idents=c("21","26"), features=rev(c("Gad1","Slc17a6","Slc18a2","Frmpd1","Lrp4","Kiss1")), cols=c("dark green", "dark blue"), ncol=1, pt.size=0)
#dev.off()
#pdf(file=paste0(output,"EstrusClustMarkers.pdf"), height=20)
#VlnPlot(mySeurat, idents=c("13"), features=rev(c("Gad1","Slc17a6","Isl1","Ebf3","Glra1","Lmo7","Abtb2","Cckar","Myo1h")), cols=c("brown"), ncol=1, pt.size=0)
#dev.off()
#a <- AverageExpression(mySeurat, features=c("Esr1", "Ar","Pgr","Esr2"))
#write.csv(a, file=paste0(output,"Esr1perclust.csv"))
typemarkers <- read.table("DotPlotMarkerLists/POA_Reordered.txt")
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
mySeurat <- ScaleData(mySeurat)
mylevels=c(0,2,7,8,9,15,16,18,22,28,30,36,37,38,1,3,4,5,6,10,11,12,13,14,17,19,20,21,23,24,25,26,27,29,31,32,33,34,35)
mylevels <- rev(mylevels)
#Idents(mySeurat) <- factor(Idents(mySeurat), levels = mylevels)
Pct <- DotPlot(mySeurat, features=typemarkers)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()
pdf(file=paste0(output,"Cyp19a1_Vln.pdf"), width=40)
VlnPlot(mySeurat, features=c("Cyp19a1"), pt.size=0)
dev.off()
pdf(file=paste0(output,"Cyp19a1_ftr.pdf"))
FeaturePlot(mySeurat, features=c("Cyp19a1"))
dev.off()
pdf(file=paste0(output,"brainstormvln.pdf"), width=40, height=60)
VlnPlot(mySeurat, features=c("Esr1","Ar","Gad1","Slc17a6","Slc18a2","Sst","Cartpt","Gal","Cck","Esr2","Tac1","Tac2","Th","Maob","Kiss1","Pappa","Col25a1","Npy","Fbn2"))
dev.off()
data <- Pct$data
head(data)
#mySeurat@active.ident <- factor(mySeurat@active.ident, levels=c("1","3","4","5","6","10","11","12","13","14","17","19","20","21","23","24","25","26","27","29","31","32","33","34","35","0","2","7","8","9","15","16","18","22","30","36","37","38"))



pdf(file=paste0(output,"_rawdotplotscaled.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()

subset <- subset(x=mySeurat, subset=Hormone=="Primed", invert=TRUE)
subset2 <- subset(x=mySeurat, subset=Hormone=="Intact", invert=TRUE)
subset3 <- subset(x=mySeurat, subset=Hormone=="Unprimed", invert=TRUE)
pdf(file=paste0(output,"Pdzrn4_Vln.pdf"), width=40)
VlnPlot(subset, features=c("Pdzrn4"), pt.size=0, idents=c("7","9","12","6","1","19","11","21","18","15","25","13","26","0","17","14","3","5","22","4","2"), split.by="Hormone", cols=c("light blue", "light green"))
dev.off()
pdf(file=paste0(output,"Pdzrncluster.pdf"))
DimPlot(mySeurat, cols=c("0"="orange", "1"="orange", "2"="orange", "3"="orange", "4"="orange", "5"="orange", "6"="orange", 
	"7"="orange", "8"="light blue", "9"="orange", "10"="light blue", "11"="orange", "12"="orange", "13"="orange", 
	"14"="orange", "15"="orange", "16"="light blue", "17"="orange", "18"="orange", "19"="orange", "20"="light blue", 
	"21"="orange", "22"="orange", "23"="light blue", "24"="light blue", "25"="orange", "26"="orange", "27"="light blue", 
	"28"="light blue", "29"="light blue", "30"="light blue", "31"="light blue", "32"="light blue", "33"="light blue", "34"="light blue","35"="light blue","36"="light blue","37"="light blue","38"="light blue"), label=FALSE, reduction="umap")
dev.off() 
pdf(file=paste0(output,"Pdzrn4_Dot.pdf"), width=40)
DotPlot(subset, features=c("Pdzrn4"), idents=c("7","9","12","6","1","19","11","21","18","15","25","13","26","0","17","14","3","5","22","4","2"), split.by="Hormone", cols=c("light blue", "light green"))
dev.off()
pdf(file=paste0(output,"Pappa_Vln.pdf"), width=20)
VlnPlot(subset, features=c("Pappa"), pt.size=0, idents=c("12","1","21","11","25","15"), split.by="Hormone", cols=c("light blue", "light green"))
dev.off()
pdf(file=paste0(output,"Mob3b_Vln.pdf"), width=20)
VlnPlot(subset, features=c("Mob3b"), pt.size=0, idents=c("1"), split.by="Hormone", cols=c("light blue", "light green"))
dev.off()
pdf(file=paste0(output,"Mob3bcluster.pdf"))
DimPlot(mySeurat, cols=c("0"="light blue", "1"="orange", "2"="light blue", "3"="light blue", "4"="light blue", "5"="light blue", "6"="light blue", 
	"7"="light blue", "8"="light blue", "9"="light blue", "10"="light blue", "11"="light blue", "12"="light blue", "13"="light blue", 
	"14"="light blue", "15"="light blue", "16"="light blue", "17"="light blue", "18"="light blue", "19"="light blue", "20"="light blue", 
	"21"="light blue", "22"="light blue", "23"="light blue", "24"="light blue", "25"="light blue", "26"="light blue", "27"="light blue", 
	"28"="light blue", "29"="light blue", "30"="light blue", "31"="light blue", "32"="light blue", "33"="light blue", "34"="light blue","35"="light blue","36"="light blue","37"="light blue","38"="light blue"), label=FALSE, reduction="umap")
dev.off() 
pdf(file=paste0(output,"Pappa_Dot.pdf"), width=20)
DotPlot(subset, features=c("Pappa"), idents=c("12","1","21","11","25","15"), split.by="Hormone", cols=c("light blue", "light green"))
dev.off()
#pdf(file=paste0(output,"Cyp19a1_Vln.pdf"), width=40)
#VlnPlot(subset, features=c("Cyp19a1"), pt.size=0, split.by="Hormone", cols=c("light blue", "light green"))
#dev.off()
pdf(file=paste0(output,"Slc17a8_Vln.pdf"))
VlnPlot(subset2, features=c("Slc17a8"), pt.size=0, idents=c("26"), split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Slc17a8cluster.pdf"))
DimPlot(mySeurat, cols=c("0"="light blue", "1"="light blue", "2"="light blue", "3"="light blue", "4"="light blue", "5"="light blue", "6"="light blue", 
	"7"="light blue", "8"="light blue", "9"="light blue", "10"="light blue", "11"="light blue", "12"="light blue", "13"="light blue", 
	"14"="light blue", "15"="light blue", "16"="light blue", "17"="light blue", "18"="light blue", "19"="light blue", "20"="light blue", 
	"21"="light blue", "22"="light blue", "23"="light blue", "24"="light blue", "25"="light blue", "26"="orange", "27"="light blue", 
	"28"="light blue", "29"="light blue", "30"="light blue", "31"="light blue", "32"="light blue", "33"="light blue", "34"="light blue","35"="light blue","36"="light blue","37"="light blue","38"="light blue"), label=FALSE, reduction="umap")
dev.off() 
pdf(file=paste0(output,"Slc17a8_Dot.pdf"))
DotPlot(subset2, features=c("Slc17a8"), idents=c("26"), split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Col25a1_Vln.pdf"), width=40)
VlnPlot(subset2, features=c("Col25a1"), pt.size=0, split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Per2_Vln.pdf"), width=40)
VlnPlot(subset2, features=c("Per2"), pt.size=0, split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Nexmif_Vln.pdf"))
VlnPlot(subset3, features=c("Nexmif"), pt.size=0, split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Acvr1c_Vln.pdf"))
VlnPlot(subset3, features=c("Acvr1c"), idents=c("11"), pt.size=0, split.by="Hormone", cols=c("blue", "pink"))
dev.off()
pdf(file=paste0(output,"Kiss1_Vln.pdf"))
VlnPlot(subset3, features=c("Kiss1"), idents=c("26"), pt.size=0, split.by="Hormone", cols=c("blue", "pink"))
dev.off()
pdf(file=paste0(output,"Kiss1cluster.pdf"))
DimPlot(mySeurat, cols=c("0"="light blue", "1"="light blue", "2"="light blue", "3"="light blue", "4"="light blue", "5"="light blue", "6"="light blue", 
	"7"="light blue", "8"="light blue", "9"="light blue", "10"="light blue", "11"="light blue", "12"="light blue", "13"="light blue", 
	"14"="light blue", "15"="light blue", "16"="light blue", "17"="light blue", "18"="light blue", "19"="light blue", "20"="light blue", 
	"21"="light blue", "22"="light blue", "23"="light blue", "24"="light blue", "25"="light blue", "26"="orange", "27"="light blue", 
	"28"="light blue", "29"="light blue", "30"="light blue", "31"="light blue", "32"="light blue", "33"="light blue", "34"="light blue","35"="light blue","36"="light blue","37"="light blue","38"="light blue"), label=FALSE, reduction="umap")
dev.off() 
pdf(file=paste0(output,"Fgfr1_Vln.pdf"), width=40)
VlnPlot(subset2, features=c("Fgfr1"), pt.size=0, split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Ybx1_Vln.pdf"), width=40)
VlnPlot(subset3, features=c("Ybx1"), pt.size=0, idents = c("3","17","25","12","15","2","0","14","6"),   split.by="Hormone", cols=c("blue", "pink"))
dev.off()
pdf(file=paste0(output,"Pid1_Vln.pdf"), width=20)
VlnPlot(subset3, features=c("Pid1"), pt.size=0, idents = c("12","11","1","22","25","9","3"),   split.by="Hormone", cols=c("light blue", "pink"))
dev.off()
pdf(file=paste0(output,"Pid1cluster.pdf"))
DimPlot(mySeurat, cols=c("0"="light blue", "1"="orange", "2"="light blue", "3"="orange", "4"="light blue", "5"="light blue", "6"="light blue", 
	"7"="light blue", "8"="light blue", "9"="orange", "10"="light blue", "11"="orange", "12"="orange", "13"="light blue", 
	"14"="light blue", "15"="light blue", "16"="light blue", "17"="light blue", "18"="light blue", "19"="light blue", "20"="light blue", 
	"21"="light blue", "22"="orange", "23"="light blue", "24"="light blue", "25"="orange", "26"="light blue", "27"="light blue", 
	"28"="light blue", "29"="light blue", "30"="light blue", "31"="light blue", "32"="light blue", "33"="light blue", "34"="light blue","35"="light blue","36"="light blue","37"="light blue","38"="light blue"), label=FALSE, reduction="umap")
dev.off() 
pdf(file=paste0(output,"Top3b_Vln.pdf"), width=40)
VlnPlot(subset2, features=c("Top3b"), pt.size=0, idents=c("2","0","21","1","6","9","3","15","16","4","12","5","17","26","13"), split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Top3bcluster.pdf"))
DimPlot(mySeurat, cols=c("0"="orange", "1"="orange", "2"="orange", "3"="orange", "4"="orange", "5"="orange", "6"="orange", 
	"7"="light blue", "8"="light blue", "9"="orange", "10"="light blue", "11"="light blue", "12"="orange", "13"="orange", 
	"14"="light blue", "15"="orange", "16"="orange", "17"="orange", "18"="light blue", "19"="light blue", "20"="light blue", 
	"21"="orange", "22"="light blue", "23"="light blue", "24"="light blue", "25"="light blue", "26"="orange", "27"="light blue", 
	"28"="light blue", "29"="light blue", "30"="light blue", "31"="light blue", "32"="light blue", "33"="light blue", "34"="light blue","35"="light blue","36"="light blue","37"="light blue","38"="light blue"), label=FALSE, reduction="umap")
dev.off() 
pdf(file=paste0(output,"Top3b_Dot.pdf"), width=40)
DotPlot(subset2, features=c("Top3b"), idents=c("2","0","21","1","6","9","3","15","16","4","12","5","17","26","13"), split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Col25a1Dotplot.pdf"))
DotPlot(mySeurat, features=c("Col25a1"), split.by="Hormone", dot.min=0.2, cols=c("blue", "red", "green"))
dev.off()
pdf(file=paste0(output,"PappaDotplot.pdf"))
DotPlot(mySeurat, features=c("Pappa"), split.by="Hormone", dot.min=0.2, cols=c("blue", "red", "green"))
dev.off()
pdf(file=paste0(output, "top3bfeatureplot.pdf"))
FeaturePlot(mySeurat, features=c("Top3b"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.25)
dev.off()
genelist <- read.table("Genelists/POA_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
MvP <- read.table("Genelists/POA_MvP_1.5cutoff.txt", header=FALSE)
MvU <- read.table("Genelists/POA_MvU_1.5cutoff.txt", header=FALSE)
sDEGs <- rbind(MvP, MvU)
sDEGs <- unlist(sDEGs)
sDEGs <- unique(sDEGs)
filtered.sDEGs <- intersect(sDEGs, genes.10x)
all <- rbind(genelist, MvP, MvU)
all <- unlist(all)
all <- unique(all)
filtered.all <- intersect(all, genes.10x)
Avg <- AverageExpression(mySeurat, assays=c("RNA"), features=filtered.genelist)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/PvUeDEGclusteravgs.csv")

pdf(file=paste0(output, "top3bfeatureplot.pdf"))
FeaturePlot(mySeurat, features=c("Top3b"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Pappafeatureplot.pdf"))
FeaturePlot(mySeurat, features=c("Pappa"),  order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Pdzrn4featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Pdzrn4"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Kiss1featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Kiss1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Mob3bfeatureplot.pdf"))
FeaturePlot(mySeurat, features=c("Mob3b"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Ybx1featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Ybx1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Acvr1cfeatureplot.pdf"))
FeaturePlot(mySeurat, features=c("Acvr1c"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Nexmiffeatureplot.pdf"))
FeaturePlot(mySeurat, features=c("Nexmif"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Pid1featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Pid1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "St3gal1featureplot.pdf"))
FeaturePlot(mySeurat, features=c("St3gal1"),  order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Slc17a8.pdf"))
FeaturePlot(mySeurat, features=c("Slc17a8"),  order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Col25a1featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Col25a1"), order=TRUE, cols=c("light blue", "red"))
dev.off()

for (i in 0:38){
ident <- i
markers <- FindMarkers(mySeurat, assay="RNA", ident.1=ident, features=filtered.all, only.pos=TRUE, min.pct=0.25, logfc.threshold=0, test.use="MAST")
write.csv(markers, file=paste0(regoutput,i,"DEGclustermarkers.csv"))
}


AvgsDEGs <- AverageExpression(mySeurat, assays=c("RNA"), features=filtered.sDEGs)
AvgsDEGs
write.csv(AvgsDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/AllsDEGclusteravgs.csv")

AvgallDEGs <- AverageExpression(mySeurat, assays=c("RNA"), features=filtered.all)
AvgallDEGs
write.csv(AvgallDEGs, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/AllDEGclusteravgs.csv")
props <- prop.table(table(Idents(mySeurat), mySeurat$Hormone), margin=2)
write.csv(props, file=paste0(output, "propperclust.csv"))
counts <- table(Idents(mySeurat), mySeurat$Hormone)
write.csv(counts, file=paste0(output,"countsperclust.csv"))
sex.markers <- FindAllMarkers(mySeurat, features=filtered.all, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.1, test.use="MAST", return.thresh=1e-4)
write.csv(sex.markers, file=paste0(output,"Sexmarkerstats.csv"))



mySeurat2 <- BuildClusterTree(mySeurat, reorder=TRUE, dims=1:30)
ScaledSeurat <- ScaleData(mySeurat2, features=rownames(mySeurat), do.center=TRUE, do.scale=TRUE)
pdf(paste0(output,"_SexMarkerheatmapreordered.pdf"))
top10 <- sex.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(ScaledSeurat, features=top10$gene, raster=FALSE) 
dev.off()
top10genes <- top10$gene
top10genes <- unique(top10genes)
Avg <- AverageExpression(mySeurat2, assays=c("RNA"), features=top10genes)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/POA/AllsigSigDEGclusteravgs.csv")
singlesig <- read.table("/scratch/users/knoedler/Heatmaps/Seurat/MetaDEGs/POA_SingleDimorphic.txt", header=FALSE)
singlesig <- unlist(singlesig)
singlesig.filtered <- setdiff(filtered.all, singlesig)

Avg2 <- AverageExpression(mySeurat, assays=c("RNA"), features=singlesig.filtered)
write.csv(Avg2, file="/scratch/users/knoedler/Heatmaps/Seurat/MetaDEGs/POAMetaDEGs.csv")

library(pheatmap)
library(viridis)
data <- read.csv("/scratch/users/knoedler/Heatmaps/Seurat/POA/AllsigSigDEGclusteravgs.csv", header=TRUE, row.names=1)
pdf(file=paste0(output,"clusteredheatmap.pdf"))
pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE, color=magma(2000))
dev.off()
