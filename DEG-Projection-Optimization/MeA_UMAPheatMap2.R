library(Seurat)
library(ggplot2)
output <- "Seurat/MeA_IndependentAnalysis/MeA_Paper_UMAPheatMap_Marchredo2"

mySeurat <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt.rds")
pdf(file=paste0(output, "Slc17a6_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Slc17a6"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Slc32a1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Slc32a1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Gad1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Gad1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Cyp19a1_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Cyp19a1"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "Nfia_featureplot.pdf"))
FeaturePlot(mySeurat, features=c("Nfia"), order=TRUE, cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output, "UMAPnobla.pdf"))
DimPlot(mySeurat, label=FALSE)
dev.off()
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
MvPSDEGs <- read.csv("Seurat/MeA_IndependentAnalysis/MeA_Merged_Clusterfuck_Marchredo_0.1logcounts_MvPSDEGS_Per_Cluster_pct.csv", header=TRUE)
MvPSDEGs <- data.frame(MvPSDEGs)
MvUSDEGs <- read.csv("Seurat/MeA_IndependentAnalysis/MeA_Merged_Clusterfuck_Marchredo_0.1logcounts_MvUSDEGS_Per_Cluster_pct.csv", header=TRUE)
MvUSDEGs <- data.frame(MvUSDEGs)
PvUEDGs <- read.csv("Seurat/MeA_IndependentAnalysis/MeA_Merged_Clusterfuck_Marchredo_0.1logcounts_PvUSDEGS_Per_Cluster_pct.csv", header=TRUE)
PvUEDGs <- data.frame(PvUEDGs)
clusters$MvP <- MvPSDEGs$MvP[match(clusters$clusters,MvPSDEGs$clusters)]
clustersMvP <- clusters["MvP"]
clusters$MvU <- MvUSDEGs$MvU[match(clusters$clusters,MvUSDEGs$clusters)]
clustersMvU <- clusters["MvU"]
clusters$PvU <- PvUEDGs$PvU[match(clusters$clusters,PvUEDGs$clusters)]
clustersPvU <- clusters["PvU"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersMvP, col.name="MvP")
mySeurat <- AddMetaData(mySeurat, metadata=clustersMvU, col.name="MvU")
mySeurat <- AddMetaData(mySeurat, metadata=clustersPvU, col.name="PvU")
pdf(file=paste0(output,"_MvPSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("MvP"), cols=c("white","black"))
dev.off()
pdf(file=paste0(output,"_MvUSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("MvU"), cols=c("white","black"))
dev.off()
pdf(file=paste0(output,"_PvUSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("PvU"), cols=c("white","black"))
dev.off()


pMvP <- FeaturePlot(mySeurat, features=c("MvP"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
pMvU <- FeaturePlot(mySeurat, features=c("MvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
pPvU <- FeaturePlot(mySeurat, features=c("PvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)

pMvPdat <- pMvP$data
pMvUdat <- pMvU$data
pPvUdat <- pPvU$data

pdf(file=paste0(output,"MvPUMAPheatmapnew.pdf"))
ggplot(pMvPdat, aes(x=UMAP_1, y=UMAP_2, color=MvP)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) + cowplot::theme_cowplot()
dev.off()

pdf(file=paste0(output,"MvUUMAPheatmapnew.pdf"))
ggplot(pMvUdat, aes(x=UMAP_1, y=UMAP_2, color=MvU)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) + cowplot::theme_cowplot()
dev.off()

pdf(file=paste0(output,"PvUUMAPheatmapnew.pdf"))
ggplot(pPvUdat, aes(x=UMAP_1, y=UMAP_2, color=PvU)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) + cowplot::theme_cowplot()
dev.off()