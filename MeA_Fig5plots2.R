#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/MeA_IndependentAnalysis/MeA_ngf_Fig5"

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
pdf(file=paste0(output,"NewUMAPsize5.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=5)
dev.off()
pdf(file=paste0(output,"NewUMAPsize6.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=6)
dev.off()
pdf(file=paste0(output,"NewUMAPsize7.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=7)
dev.off()
pdf(file=paste0(output,"NewUMAPsize8.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=8)
dev.off()
pdf(file=paste0(output,"UMAP_Hormonelabelshuffle.pdf"))
DimPlot(mySeurat, reduction="umap", group.by="Hormone", shuffle=TRUE, seed=1, label=FALSE, pt.size=0.4, cols=c("blue", "red", "green"))
dev.off()
pdf(file=paste0(output,"brainstormvln.pdf"), width=40, height=80)
VlnPlot(mySeurat, features=c("Esr1","Ar","Gad1","Slc17a6","Slc18a2","Sst","Cartpt","Gal","Cck","Esr2","Tac1","Tac2","Th","Maob","Kiss1","Pappa","Col25a1","Npy","Fbn2"), ncol=1, pt.size=0)
dev.off()
pdf(file=paste0(output,"tonsomarkers_ftr.pdf"), width=50, height=50)
FeaturePlot(mySeurat, features=c("Gad1","Slc32a1","Gad2","Slc17a6","Slc1a3","Slc29a4","Nfia","Nfix","Nfib","Slc17a7","Esr1"), cols=c("light blue","red"))
dev.off()
pdf(file=paste0(output,"neurlin.pdf"), width=40, height=100)
VlnPlot(mySeurat, features=c("Syn1","Syn2","Syn3","Rbfox3","Nfia","Nfib","Nfix","Eno1","Syp","Map2","Slc1a3","Slc29a4","Slc17a7","Reln","Nos1","Npy","Nr2f2","Car12","Pparg","Met","Aqp4","Nes","Gfap","Olig2","Rax","Sox2","Sox10","Vim","Hepacam","Plp1","Slc4a4"), ncol=1, pt.size=0)
dev.off()
pdf(file=paste0(output,"Gad1_ftr.pdf"))
FeaturePlot(mySeurat, features=c("Gad1"), order=TRUE, cols=c("light blue","red"))
dev.off()
pdf(file=paste0(output,"Slc17a6_ftr.pdf"))
FeaturePlot(mySeurat, features=c("Slc17a6"), order=TRUE, cols=c("light blue","red"))
dev.off()
pdf(file=paste0(output,"Esr1_ftr.pdf"))
FeaturePlot(mySeurat, features=c("Esr1"), order=TRUE, cols=c("light blue","red"))
dev.off()
ngf <- subset(mySeurat, idents=c("4","9","25"))
ngf.markers <- FindAllMarkers(ngf, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.30, test.use="MAST",  min.diff.pct = 0.2)
write.csv(ngf.markers, file=paste0(output,"ngfmarkers.csv"))
ngf <- ScaleData(ngf, features=rownames(ngf), do.center=TRUE, do.scale=TRUE)
mySeurat <- ScaleData(mySeurat, features=rownames(ngf), do.center=TRUE, do.scale=TRUE)
pdf(paste0(output,"_Ngftop10Markerheatmap.pdf"))
top10 <- ngf.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(ngf, features=top10$gene, raster=FALSE) 
dev.off()

pdf(file=paste0(output,"SFARIheatmap.pdf"), width=20)
DoHeatmap(mySeurat, features=c("Auts2","Tcf4","Sox5","Adnp","Nr3c2","Cntnap2"))
dev.off()
pdf(file=paste0(output,"NGF_SFARIheatmap.pdf"), width=20)
DoHeatmap(ngf, features=c("Auts2","Tcf4","Sox5","Adnp","Nr3c2","Cntnap2"))
dev.off()

top10genes <- top10$gene
top10genes <- unique(top10genes)
Avg <- AverageExpression(ngf, assays=c("RNA"), features=top10genes)
Avg
write.csv(Avg, file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/NgfMarker_top10clusteravgs.csv")
pdf(paste0(output,"_Ngftop5Markerheatmap.pdf"))
top5 <- ngf.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
DoHeatmap(ngf, features=top5$gene, raster=FALSE) 
dev.off()

top5genes <- top5$gene
top5genes <- unique(top5genes)
Avg5 <- AverageExpression(ngf, assays=c("RNA"), features=top5genes)
Avg5
write.csv(Avg5, file="/scratch/users/knoedler/Heatmaps/Seurat/MeA/NgfMarker_top5clusteravgs.csv")

pdf(file=paste0(output,"_Col25a1_hormonesplit.pdf"), width=15)
VlnPlot(ngf, features=c("Col25a1"), split.by="Hormone", cols=c("light blue", "light pink", "light green"), pt.size=0)
dev.off()
pdf(file=paste0(output,"_Sox5_hormonesplit.pdf"), width=15)
VlnPlot(ngf, features=c("Sox5"), split.by="Hormone", cols=c("light blue", "light pink", "light green"), pt.size=0)
dev.off()
pdf(file=paste0(output,"_Chrna7_hormonesplit.pdf"), width=15)
VlnPlot(ngf, features=c("Chrna7"), split.by="Hormone", cols=c("light blue", "light pink", "light green"), pt.size=0)
dev.off()

pdf(file=paste0(output,"9_Col25a1_hormonesplit.pdf"))
VlnPlot(ngf, idents=c("9"), features=c("Col25a1"), split.by="Hormone", cols=c("light blue", "light pink", "light green"), pt.size=0)
dev.off()
pdf(file=paste0(output,"9_Sox5_hormonesplit.pdf"))
VlnPlot(ngf,  idents=c("9"), features=c("Sox5"), split.by="Hormone", cols=c("light blue", "light pink", "light green"), pt.size=0)
dev.off()
pdf(file=paste0(output,"9_Chrna7_hormonesplit.pdf"))
VlnPlot(ngf,  idents=c("9"), features=c("Chrna7"), split.by="Hormone", cols=c("light blue", "light pink", "light green"), pt.size=0)
dev.off()

mySeurat <- RenameIdents(mySeurat, "0"="excit")
mySeurat <- RenameIdents(mySeurat, "3"="excit")
mySeurat <- RenameIdents(mySeurat, "8"="excit")
mySeurat <- RenameIdents(mySeurat, "12"="excit")
mySeurat <- RenameIdents(mySeurat, "13"="excit")
mySeurat <- RenameIdents(mySeurat, "14"="excit")
mySeurat <- RenameIdents(mySeurat, "16"="excit")
mySeurat <- RenameIdents(mySeurat, "17"="excit")
mySeurat <- RenameIdents(mySeurat, "20"="excit")
mySeurat <- RenameIdents(mySeurat, "21"="excit")
mySeurat <- RenameIdents(mySeurat, "22"="excit")
mySeurat <- RenameIdents(mySeurat, "28"="excit")
mySeurat <- RenameIdents(mySeurat, "30"="excit")
mySeurat <- RenameIdents(mySeurat, "31"="excit")
mySeurat <- RenameIdents(mySeurat, "32"="excit")
mySeurat <- RenameIdents(mySeurat, "33"="excit")
#excit <- c("0","3","8","12","13","14","16","17","20","21","22","28","30","31","32","33")
#excit
#excit <- unlist(excit)
#excit
#inhib <- c("1","2","5","6","7","10","16","15","18","19","23","24","26","27","28","29","34")
mySeurat <- RenameIdents(mySeurat, "1"="inhib")
mySeurat <- RenameIdents(mySeurat, "2"="inhib")
mySeurat <- RenameIdents(mySeurat, "5"="inhib")
mySeurat <- RenameIdents(mySeurat, "6"="inhib")
mySeurat <- RenameIdents(mySeurat, "7"="inhib")
mySeurat <- RenameIdents(mySeurat, "10"="inhib")
mySeurat <- RenameIdents(mySeurat, "11"="inhib")
mySeurat <- RenameIdents(mySeurat, "15"="inhib")
mySeurat <- RenameIdents(mySeurat, "18"="inhib")
mySeurat <- RenameIdents(mySeurat, "19"="inhib")
mySeurat <- RenameIdents(mySeurat, "23"="inhib")
mySeurat <- RenameIdents(mySeurat, "24"="inhib")
mySeurat <- RenameIdents(mySeurat, "26"="inhib")
mySeurat <- RenameIdents(mySeurat, "27"="inhib")
mySeurat <- RenameIdents(mySeurat, "29"="inhib")
mySeurat <- RenameIdents(mySeurat, "34"="inhib")
mySeurat <- RenameIdents(mySeurat, "4"="NGF")
mySeurat <- RenameIdents(mySeurat, "9"="NGF")
mySeurat <- RenameIdents(mySeurat, "25"="NGF")
#inhib <- unlist(inhib)
#for (g in inhib){
#old <- g
#mySeurat=RenameIdents(mySeurat, old.ident.name = old, new.ident.name='inhib')
#}
#for(g in excit){
#old <- g
#mySeurat <- RenameIdents(mySeurat, old.ident.name = old, new.ident.name='excit')
#}

pdf(file=paste0(output,"MarkerPlot.pdf"))
DotPlot(mySeurat, features=c("Slc17a7","Slc1a3","Slc4a4","Nfia","Nfib","Nfix"))
dev.off()
p <- DotPlot(mySeurat, features=c("Slc17a7","Slc1a3","Slc4a4","Nfia","Nfib","Nfix"))
data <- p$data
pdf(file=paste0(output,"_ngfplotscaled.pdf"))
data %>% 
  filter(pct.exp > 10) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + 
  scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())   + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()