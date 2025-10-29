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

output <- "Seurat/BNST_IndependentAnalysis/BNST_Figs5heatmapsfinal_markers_refilteredfinal"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/BNST_MvP_1.5cutoff.txt")
eDEG <- read.table("Genelists/BNST_PvU_1.5cutoff.txt")
sDEG2 <- read.table("Genelists/BNST_MvU_1.5cutoff.txt")
unlist(eDEG)
unlist(sDEG2)
unlist(Dimorphic)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered.m <- as.matrix(eDEG.filtered)

sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered.m <- as.matrix(sDEG2.filtered)

Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered.m <- as.matrix(Dimorphic.filtered)
mySeurat <- ScaleData(mySeurat)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
mylevels <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)
mylevels <- rev(mylevels)
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)
specDEGs <- c("Rbm20","Plekhg1","Prok2","Cck","Phf21b","Slco2a1","Tac1","Asic4","Kank1","Sst","Abca1","Amot","Fam107b","Prkcd","Sertm1","Tns1","Sytl4","Npr3","Plscr4")
Pct <- DotPlot(mySeurat, features=specDEGs)

data <- Pct$data

pdf(file=paste0(output,"_specDEGdotplot_scaledatabluered.pdf"))
data %>% 
  ggplot(aes(x=features.plot, y=factor(id, levels=mylevels), color=avg.exp.scaled, size=pct.exp)) + 
  geom_point() + scale_colour_gradientn(colours = c("blue","red"), limits=c(-1.5, 2.5))+ scale_size(limits=c(25,100), breaks=seq(25,100,25)) +
  cowplot::theme_cowplot() + theme(axis.title.x=element_text(color="black", size=20)) + labs(x="") +
  theme(axis.text.x = element_text(angle = 45)) +
  ylab('') + scale_x_discrete(position = "top") +
  theme(axis.ticks = element_blank())
dev.off()

AllsDegs <- read.table("topGO/Total_ByRegion/BNST_Genesonly.txt")
AllsDegs <- unlist(AllsDegs)
AllsDegs.m <- as.matrix(AllsDegs)
AllsDegs.filtered <- intersect(AllsDegs, genes.10x)
AllsDegs.filtered.m <- as.matrix(AllsDegs.filtered)
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
AvgFinal <- AverageExpression(mySeurat, features=c("Rbm20","Plekhg1","Prok2","Cck","Phf21b","Slco2a1","Tac1","Asic4","Kank1","Sst","Abca1","Amot","Fam107b","Prkcd","Sertm1","Tns1","Sytl4","Npr3","Plscr4"))
AvgFinal2 <- AvgFinal$RNA
ScaledFinal <- scale(t(AvgFinal2))
write.csv(ScaledFinal, file=paste0(output,"BNST_Sexmarkers_Allclusters.csv"))
pdf(file=paste0(output,"Allclusters_UniqueSexmarkers.pdf"))
pheatmap(ScaledFinal, cluster_rows=FALSE, cluster_cols=FALSE, color=magma(2000))
dev.off()
#Markers <- FindAllMarkers(mySeurat, assay="RNA", only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.20)
#write.csv(Markers, file=paste0(output,"allMarkerswtfwilcox.csv"))
All.markers <- FindAllMarkers(mySeurat, assay="RNA", features=AllsDegs.filtered, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.223, test.use="MAST")
write.csv(All.markers, file=paste0(output,"_RNA.csv"))
#saveRDS(Markers, file=paste0(output,"precomputed.rds"))
Sig.Markers <- All.markers[All.markers$p_val_adj < 0.05,]
Sig.Markers <- Sig.Markers[Sig.Markers$pct.2 < 0.25,]
head(Sig.Markers)
siggenes <- Sig.Markers$gene
head(siggenes)
uniquesig <- names(table(siggenes)[table(siggenes)==1])
uniquesig
uniquesig <- unlist(uniquesig)
uniquesig.m <- as.matrix(uniquesig)
uniqueclustssig <- Sig.Markers[uniquesig.m,]
uniqueclustssig
MvP <- uniqueclustssig[Dimorphic.filtered.m,]
MvU <- uniqueclustssig[sDEG2.filtered.m,]
PvU <- uniqueclustssig[eDEG.filtered.m,]
MarkedMvP <- MvP$cluster
MarkedMvP <- unique(MarkedMvP)
MarkedMvP <- na.omit(MarkedMvP)
MarkedMvU <- MvU$cluster
MarkedMvU <- unique(MarkedMvU)
MarkedPvU <- PvU$cluster
MarkedPvU <- unique(PvU)
write.csv(uniqueclustssig, file=paste0(output,"MvP_uniquemarkerswithclusters.csv"))
Avg <- AverageExpression(mySeurat, assays=c("RNA"), features=uniquesig)
Avg
write.csv(Avg, file=paste0(output,"allsigUniquemarkers_AvgExpression.csv"))
names(Avg) = gsub(pattern="RNA.", replacement="", x=names(Avg))
Avg
Avg2 <- Avg$RNA
Avg2
MarkedmCTssig <- uniqueclustssig$cluster
MarkedmCTssig  <- unique(MarkedmCTssig)
MarkedAvg <- Avg2[,MarkedmCTssig]
MarkedAvg
MarkedAvg <- scale(t(MarkedAvg))
pdf(file=paste0(output,"allSigmarkerDEGs.pdf"))
pheatmap(MarkedAvg, col=magma(2000))
dev.off()
write.csv(MarkedAvg, file=paste0(output,"SexMarkersscaled.csv"))
MvPavg <- Avg2[Dimorphic.filtered.m,]
pdf(file=paste0(output, "Dotplot25.pdf"), width=20)
DotPlot(mySeurat, features=uniquesig, dot.min=0.25)
dev.off()
pdf(file=paste0(output, "Dotplot11.pdf"), width=20)
DotPlot(mySeurat, features=uniquesig, dot.min=0.10)
dev.off()
