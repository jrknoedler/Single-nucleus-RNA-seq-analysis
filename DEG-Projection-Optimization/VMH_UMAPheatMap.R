library(Seurat)
library(ggplot2)
output <- "Seurat/VMH_IndependentAnalysis/VMH_Paper_UMAPheatMap_"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_filtered3_Paperanalysis_20pcs.rds")

clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
MvPSDEGs <- read.csv("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_20pcs_0.1cutoff_1.5foldMvP_SDEGS_Per_Cluster.csv", header=TRUE)
MvPSDEGs <- data.frame(MvPSDEGs)
MvUSDEGs <- read.csv("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_20pcs_0.1cutoff_1.5foldMvU_SDEGS_Per_Cluster.csv", header=TRUE)
MvUSDEGs <- data.frame(MvUSDEGs)
PvUEDGs <- read.csv("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_20pcs_0.1cutoff_1.5foldPvU_SDEGS_Per_Cluster.csv", header=TRUE)
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
FeaturePlot(mySeurat, features=c("MvP"), cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output,"_MvUSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("MvU"), cols=c("light blue", "red"))
dev.off()
pdf(file=paste0(output,"_PvUSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("PvU"), cols=c("light blue", "red"))
dev.off()
embeddings <- Embeddings(mySeurat, reduction="umap")
embeddings <- data.frame(embeddings)
MvPq <- mySeurat[["MvP"]]
MvPq <- data.frame(MvPq)
MvUq <- mySeurat[["MvU"]]
MvUq <- data.frame(MvUq)
PvUq <- mySeurat[["PvU"]]
PvUq <- data.frame(PvUq)
UMAP1 <- merge(embeddings, MvPq)
head(UMAP1)
UMAP2 <- merge(embeddings, MvUq)
head(UMAP2)
UMAP3 <- merge(embeddings, PvUq)
head(UMAP3)
pdf(file=paste0(output,"_MvP_EditedUMAP.pdf"))
plot <- ggplot(UMAP1, aes(x=UMAP_1, y=UMAP_2, color=MvP)) + geom_point(aes(color=MvP)) + theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ theme(plot.background= element_rect(fill="white"))
plot
dev.off()
pdf(file=paste0(output,"_MvU_EditedUMAP.pdf"))
plot <- ggplot(UMAP2, aes(x=UMAP_1, y=UMAP_2, color=MvU)) + geom_point(aes(color=MvP)) + theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ theme(plot.background= element_rect(fill="white"))
plot
dev.off()
pdf(file=paste0(output,"_PvU_EditedUMAP.pdf"))
plot <- ggplot(UMAP3, aes(x=UMAP_1, y=UMAP_2, color=PvU)) + geom_point(aes(color=MvP)) + theme(panel.background = element_rect(fill = "white", color = "white", linetype="solid")) + theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank())+ theme(plot.background= element_rect(fill="white"))
plot
dev.off()