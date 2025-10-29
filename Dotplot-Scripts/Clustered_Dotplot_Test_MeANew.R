#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/MeA_Betterredofilt_lowdot"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(tidyverse)
library(ggtree)



mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")
pdf(file=paste0(output,"UMAP_Hormonelabelshuffle.pdf"))
DimPlot(mySeurat, reduction="umap", group.by="Hormone", shuffle=TRUE, seed=1, label=FALSE, pt.size=0.4, cols=c("blue", "red", "green"))
dev.off()
pdf(file=paste0(output,"brainstormvln.pdf"), width=40, height=80)
VlnPlot(mySeurat, features=c("Esr1","Ar","Gad1","Slc17a6","Slc18a2","Sst","Cartpt","Gal","Cck","Esr2","Tac1","Tac2","Th","Maob","Kiss1","Pappa","Col25a1","Npy","Fbn2"), ncol=1, pt.size=0)
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
a <- AverageExpression(mySeurat, assays=c("RNA"), features=c("Esr1", "Ar","Pgr","Esr2"))
write.csv(a, file=paste0(output,"Esr1perclust.csv"))
pdf(file=paste0(output,"Nfiavln.pdf"), width=20, height=10)
VlnPlot(mySeurat, features=c("Nfia","Nfix","Auts2","Nr3c2","Slc1a3"), pt.size=0, split.by="Hormone", ncol=1)
dev.off()
pdf(file=paste0(output,"biasmarkers.pdf"), width=10, height=20)
VlnPlot(mySeurat, idents=c("1","13","18","21","23","24","27","29"), features=rev(c("Slc1a3","Slc17a6","Adamts2","Piezo2","Ebf3","Frmpd1","Chd7","Gm45194","Glra1","Col5a1")), ncol=1, pt.size=0)
dev.off()
subset <- subset(x=mySeurat, subset=Hormone=="Primed", invert=TRUE)
subset2 <- subset(x=mySeurat, subset=Hormone=="Intact", invert=TRUE)
subset3 <- subset(x=mySeurat, subset=Hormone=="Unprimed", invert=TRUE)
pdf(file=paste0(output,"NfiaMvPvln.pdf"), width=40)
VlnPlot(subset3, features=c("Nfia"), pt.size=0, split.by="Hormone", cols=c("light blue", "pink"))
dev.off()
pdf(file=paste0(output,"Cyp19a1MvPvln.pdf"), width=40)
VlnPlot(subset3, features=c("Cyp19a1"), pt.size=0, split.by="Hormone", cols=c("light blue", "pink"))
dev.off()
pdf(file=paste0(output,"Auts2MvPvln.pdf"), width=40)
VlnPlot(subset3, features=c("Auts2"), pt.size=0, split.by="Hormone", cols=c("light blue", "pink"))
dev.off()
pdf(file=paste0(output,"NfixMvPvln.pdf"), width=40)
VlnPlot(subset3, features=c("Nfix"), pt.size=0, split.by="Hormone", cols=c("light blue", "pink"))
dev.off()
pdf(file=paste0(output,"Cntnap3MvUvln.pdf"), width=40)
VlnPlot(subset, features=c("Cntnap3"), pt.size=0, split.by="Hormone", cols=c("light blue", "light green"))
dev.off()
pdf(file=paste0(output,"Spata13MvUvln.pdf"), width=40)
VlnPlot(subset, features=c("Spata13"), pt.size=0, split.by="Hormone", cols=c("light blue", "light green"))
dev.off()
pdf(file=paste0(output,"Npas1PvUvln.pdf"), width=40)
VlnPlot(subset2, features=c("Npas2"), pt.size=0, split.by="Hormone", , cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Col25a1PvUvln.pdf"), width=40)
VlnPlot(subset2, features=c("Col25a1"), pt.size=0, split.by="Hormone", cols=c("pink", "light green"))
dev.off()
pdf(file=paste0(output,"Clusterbias.pdf"))
DimPlot(mySeurat, cols=c("0"="light gray", "1"="purple", "2"="light gray", "3"="light gray", "4"="light gray", "5"="light gray", "6"="light gray", 
	"7"="light gray", "8"="light gray", "9"="light gray", "10"="light gray", "11"="light gray", "12"="light gray", "13"="light green", 
	"14"="light gray", "15"="light gray", "16"="light gray", "17"="light gray", "18"="blue", "19"="light gray", "20"="light gray", 
	"21"="orange", "22"="light gray", "23"="purple", "24"="purple", "25"="light gray", "26"="light gray", "27"="purple", 
	"28"="light gray", "29"="purple", "30"="light gray", "31"="light gray", "32"="light gray", "33"="light gray","34"="light gray"), label=FALSE, reduction="umap")
dev.off() 
typemarkers <- read.table("DotPlotMarkerLists/MeA_Reorganized_Filtered3_CCAmerged.txt")
pdf(file=paste0(output,"plotnolab.pdf"))
DimPlot(mySeurat, label=FALSE, reduction="umap")
dev.off()
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
mySeurat <- NormalizeData(mySeurat)
mySeurat <- ScaleData(mySeurat)
Pct <- DotPlot(mySeurat, features=typemarkers)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
pdf(file=paste0(output,"_rawdotplotscaled.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()


# make data square to calculate euclidean distance
#mat <- data %>% 
#  select(-cluster) %>%  # drop unused columns to faciliate widening
#  pivot_wider(names_from = cluster, values_from = count) %>% 
#  data.frame() # make df as tibbles -> matrix annoying
#row.names(mat) <- mat$Gene  # put gene in `row`
#mat <- mat[,-1] #drop gene column as now in rows
#clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix

#pdf(file=paste0(output,"clusteredDotplot"))
#dotplot <- gene_cluster %>% 
# # mutate(Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
#  #filter(pct.exp > 25) %>% 
#  ggplot(aes(x=id, y=features.plot, color=avg.exp, size=pct.exp))  + 
#  geom_point() + 
#  cowplot::theme_cowplot() + 
#  theme(axis.line  = element_blank()) +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
##  ylab('') +
#  theme(axis.ticks = element_blank()) +
#  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')
#dev.off()
genelist <- read.table("Genelists/MeA_PvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
MvP <- read.table("Genelists/MeA_MvP_1.5cutoff.txt", header=FALSE)
MvU <- read.table("Genelists/MeA_MvU_1.5cutoff.txt", header=FALSE)
sDEGs <- rbind(MvP, MvU)
sDEGs <- unlist(sDEGs)
sDEGs <- unique(sDEGs)
filtered.sDEGs <- intersect(sDEGs, genes.10x)
all <- rbind(genelist, MvP, MvU)
all <- unlist(all)
all <- unique(all)
filtered.all <- intersect(all, genes.10x)
mySeurat <- NormalizeData(mySeurat)
props <- prop.table(table(Idents(mySeurat), mySeurat$Hormone), margin=2)
write.csv(props, file=paste0(output, "propperclust.csv"))
counts <- table(Idents(mySeurat), mySeurat$Hormone)
write.csv(counts, file=paste0(output,"countsperclust.csv"))
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
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/AllsigDEGclusteravgs.csv")

singlesig <- read.table("/scratch/users/knoedler/Heatmaps/Seurat/MetaDEGs/MeA_SingleDimorphic.txt", header=FALSE)
singlesig <- unlist(singlesig)
singlesig.filtered <- setdiff(filtered.all, singlesig)

Avg2 <- AverageExpression(mySeurat, assays=c("RNA"), features=singlesig.filtered)
write.csv(Avg2, file="/scratch/users/knoedler/Heatmaps/Seurat/MetaDEGs/MeAMetaDEGs.csv")

library(pheatmap)
data <- read.csv("/scratch/users/knoedler/Heatmaps/Seurat/MeA/AllsigDEGclusteravgs.csv", header=TRUE, row.names=1)
pdf(file=paste0(output,"clusteredheatmap.pdf"))
pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE, color=magma(2000))
dev.off()



pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.25)
dev.off()
genelist <- read.table("Genelists/MeA_MvU_1.5cutoff.txt", header=FALSE)
genelist <- unlist(genelist)
dim(x=mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
filtered.genelist <- intersect(genelist, genes.10x)
filtered.genelist
mySeurat <- NormalizeData(mySeurat)
Avg <- AverageExpression(mySeurat, assays=c("RNA"), features=filtered.genelist)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/MvUsDEGclusteravgs.csv")
