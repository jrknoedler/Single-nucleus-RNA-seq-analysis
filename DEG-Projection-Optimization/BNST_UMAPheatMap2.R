library(Seurat)
library(ggplot2)
output <- "Seurat/BNST_IndependentAnalysis/BNST_Paper_UMAPheatMap_grayscale"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Esr1filtered_striatumEsr1filtered_malatexcludefiltered_sexclude_malat1regress_30pcsres1.2.rds")
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)

clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
MvPSDEGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Marchredo_0.5logcounts_MvPSDEGS_Per_Cluster_pct.csv", header=TRUE)
MvPSDEGs <- data.frame(MvPSDEGs)
MvUSDEGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Marchredo_0.5logcounts_MvUSDEGS_Per_Cluster_pct.csv", header=TRUE)
MvUSDEGs <- data.frame(MvUSDEGs)
PvUEDGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Marchredo_0.5logcounts_PvUSDEGS_Per_Cluster_pct.csv", header=TRUE)
PvUEDGs <- data.frame(PvUEDGs)

sigDEGs <- read.table("ClusterDeg_Forplots/BNST.txt", header=TRUE, sep="\t")
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

pdf(file=paste0(output,"_BNST_SigDegsperclust.pdf"), width=40, height=50)
FeaturePlot(mySeurat, features=c("MvPsig","MvUsig","PvUSig","sDEGSig","TotalDegSig","Esr1"))
dev.off()

pdf(file=paste0(output,"_MvPSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("MvP"), cols=c("white","black")) 
dev.off()
pdf(file=paste0(output,"_MvUSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("MvU"), cols=c("white","black"))
dev.off()
pdf(file=paste0(output,"_PvUSDegsperclust.pdf")) 
FeaturePlot(mySeurat, features=c("PvU"), col=c("white","black"))
dev.off()

pMvP <- FeaturePlot(mySeurat, features=c("MvP"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
pMvU <- FeaturePlot(mySeurat, features=c("MvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
pPvU <- FeaturePlot(mySeurat, features=c("PvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)

pMvPsig <- FeaturePlot(mySeurat, features=c("MvPsig"))
pMvUsig <- FeaturePlot(mySeurat, features=c("MvUsig"))
pPvUsig <- FeaturePlot(mySeurat, features=c("PvUSig"))
psDEGsig <- FeaturePlot(mySeurat, features=c("sDEGSig"))
pTotalsig <- FeaturePlot(mySeurat, features=c("TotalDegSig"))

pdf(file=paste0(output,"Esr1ftr.pdf"))
FeaturePlot(mySeurat, features=c("Esr1"), cols=c("light blue", "orange"))
dev.off()

pMvPdat <- pMvP$data
pMvUdat <- pMvU$data
pPvUdat <- pPvU$data
pMvPsigdat <- pMvPsig$data
pMvUsigdat <- pMvUsig$data
pPvUsigdat <- pPvUsig$data
psDEGsigdat <- psDEGsig$data
pTotalsigdat <- pTotalsig$data

pdf(file=paste0(output,"MvPsigUMAPheatmapnew.pdf"))
ggplot(pMvPsigdat, aes(x=UMAP_1, y=UMAP_2, color=MvPsig)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("light blue","red"), limits=c(0, 185)) + cowplot::theme_cowplot()
dev.off()
pdf(file=paste0(output,"MvUsigUMAPheatmapnew.pdf"))
ggplot(pMvUsigdat, aes(x=UMAP_1, y=UMAP_2, color=MvUsig)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("light blue","red"), limits=c(0, 185)) + cowplot::theme_cowplot()
dev.off()
pdf(file=paste0(output,"PvUsigUMAPheatmapnew.pdf"))
ggplot(pPvUsigdat, aes(x=UMAP_1, y=UMAP_2, color=PvUSig)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("gray","red"), limits=c(0, 35)) + cowplot::theme_cowplot()
dev.off()
pdf(file=paste0(output,"sDEGsigUMAPheatmapnew.pdf"))
ggplot(psDEGsigdat, aes(x=UMAP_1, y=UMAP_2, color=sDEGSig)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("gray","blue"), limits=c(0, 160)) + cowplot::theme_cowplot()
dev.off()
pdf(file=paste0(output,"TotalsigUMAPheatmapnew.pdf"))
ggplot(pTotalsigdat, aes(x=UMAP_1, y=UMAP_2, color=TotalDegSig)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("light blue","red"), limits=c(0, 185)) + cowplot::theme_cowplot()
dev.off()

pdf(file=paste0(output,"MvPUMAPheatmapnew.pdf"))
ggplot(pMvPdat, aes(x=UMAP_1, y=UMAP_2, color=MvP)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) + cowplot::theme_cowplot()
dev.off()

pdf(file=paste0(output,"MvUUMAPheatmapnew.pdf"))
ggplot(pMvUdat, aes(x=UMAP_1, y=UMAP_2, color=MvU)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) + cowplot::theme_cowplot()
dev.off()

pdf(file=paste0(output,"PvUUMAPheatmapnew.pdf"))
ggplot(pPvUdat, aes(x=UMAP_1, y=UMAP_2, color=PvU)) +geom_point(size=0.3) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) + cowplot::theme_cowplot()
dev.off()