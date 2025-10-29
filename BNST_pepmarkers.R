#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"


library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Figs6dotplot_reclust"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")

DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
typemarkers <- read.csv("Genelists/Supp6_7_dotplots.csv", header=FALSE)
typemarkers <- unlist(typemarkers)
typemarkers.filtered <- intersect(typemarkers, genes.10x)
head(typemarkers.filtered)
typemarkers.filtered <- as.matrix(typemarkers.filtered)
#Markers <- FindAllMarkers(mySeurat, assay="RNA", features=typemarkers.filtered, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.20)
#write.csv(Markers, file=paste0(output,"allMarkerswtfwilcox.csv"))
All.markers <- FindAllMarkers(mySeurat, assay="RNA", features=typemarkers.filtered, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.223, test.use="MAST")
write.csv(All.markers, file=paste0(output,"_RNA.csv"))
#saveRDS(Markers, file=paste0(output,"precomputed.rds"))
Sig.Markers <- All.markers[All.markers$p_val_adj < 0.05,]
Sig.Peptides <- Sig.Markers[typemarkers.filtered,]
Pepgenes <- Sig.Peptides$gene
Pepgenes <- unique(Pepgenes)
mylevels <- c(6,10,20,33,0,1,2,3,4,5,7,8,9,11,12,13,14,15,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32)
mylevels <- rev(mylevels)
Pct <- DotPlot(mySeurat, features=Pepgenes)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
data <- Pct$data




pdf(file=paste0(output,"_filtereddotplot_scaledatabluered.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"))+ scale_size(limits=c(10,100)) +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())  + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()


