library(Seurat)
library(ggplot2)
output <- "Seurat/POA_IndependentAnalysis/POA_Paper_UMAPheatMap_30pcs_grayscaleMarchredo"

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")

clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
MvPSDEGs <- read.csv("Seurat/POA_IndependentAnalysis/POA_Merged_Clusterfuck_Marchredo_0.1logcounts_MvPSDEGS_Per_Cluster_pct.csv", header=TRUE)
MvPSDEGs <- data.frame(MvPSDEGs)
MvUSDEGs <- read.csv("Seurat/POA_IndependentAnalysis/POA_Merged_Clusterfuck_Marchredo_0.1logcounts_MvUSDEGS_Per_Cluster_pct.csv", header=TRUE)
MvUSDEGs <- data.frame(MvUSDEGs)
PvUEDGs <- read.csv("Seurat/POA_IndependentAnalysis/POA_Merged_Clusterfuck_Marchredo_0.1logcounts_PvUSDEGS_Per_Cluster_pct.csv", header=TRUE)
PvUEDGs <- data.frame(PvUEDGs)
sigDEGs <- read.table("ClusterDeg_Forplots/POA.txt", header=TRUE, sep="\t")
sigDEGs <- data.frame(sigDEGs)
sigDEGs
sigDEGsMvP <- sigDEGs[,c(1,2)]
sigDEGsMvP
sigDEGsMvU <- sigDEGs[,c(1,3)]
sigDEGsPvU <- sigDEGs[,c(1,4)]
sigDEGsTotals <- sigDEGs[,c(1,5)]
sigDEGsTotald <- sigDEGs[,c(1,6)]

clusters$MvP <- MvPSDEGs$MvP[match(clusters$clusters,MvPSDEGs$clusters)]
clustersMvP <- clusters["MvP"]
clusters$MvU <- MvUSDEGs$MvU[match(clusters$clusters,MvUSDEGs$clusters)]
clustersMvU <- clusters["MvU"]
clusters$PvU <- PvUEDGs$PvU[match(clusters$clusters,PvUEDGs$clusters)]
clustersPvU <- clusters["PvU"]

clusters$MvPSig <- sigDEGsMvP$MvP[match(clusters$clusters,sigDEGsMvP$clusters)]
clustersMvPSig <- clusters["MvPSig"]

clusters$MvUSig <- sigDEGsMvU$MvU[match(clusters$clusters,sigDEGsMvU$clusters)]
clustersMvUSig <- clusters["MvUSig"]

clusters$PvUSig <- sigDEGsPvU$PvU[match(clusters$clusters,sigDEGsPvU$clusters)]
clustersPvUSig <- clusters["PvUSig"]

clusters$TotalsDEG <- sigDEGsTotals$TotalsDEG[match(clusters$clusters,sigDEGsTotals$clusters)]
clusterssDEGSig <- clusters["TotalsDEG"]

clusters$TotalSig <- sigDEGsTotald$TotalDEG[match(clusters$clusters,sigDEGsTotald$clusters)]
clustersDEGSig <- clusters["TotalSig"]

mySeurat <- AddMetaData(mySeurat, metadata=clustersMvP, col.name="MvP")
mySeurat <- AddMetaData(mySeurat, metadata=clustersMvU, col.name="MvU")
mySeurat <- AddMetaData(mySeurat, metadata=clustersPvU, col.name="PvU")
mySeurat <- AddMetaData(mySeurat, metadata=clustersMvPSig, col.name="MvPsig")
mySeurat <- AddMetaData(mySeurat, metadata=clustersMvUSig, col.name="MvUsig")
mySeurat <- AddMetaData(mySeurat, metadata=clustersPvUSig, col.name="PvUSig")
mySeurat <- AddMetaData(mySeurat, metadata=clusterssDEGSig, col.name="sDEGSig")
mySeurat <- AddMetaData(mySeurat, metadata=clustersDEGSig, col.name="TotalDegSig")

pdf(file=paste0(output,"_POA_SigDegsperclust.pdf"), width=40, height=50)
FeaturePlot(mySeurat, features=c("MvPsig","MvUsig","PvUSig","sDEGSig","TotalDegSig","Esr1"))
dev.off()


pdf(file=paste0(output,"_POA_SDegsperclust.pdf"), width=16, height=16)
FeaturePlot(mySeurat, features=c("MvP","MvU","PvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
dev.off()
p <- FeaturePlot(mySeurat, features=c("MvP","MvU","PvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
head(p$data)


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

p <- DimPlot(mySeurat, reduction="umap", group.by="Hormone")
data <- p$data
head(data)
write.csv(data, file=paste0(output,"forChungha.csv"))

