library(Seurat)
library(ggplot2)
output <- "Seurat/POA_IndependentAnalysis/POA_Paper_UMAPheatMap_30pcs_grayscale"

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test.rds")

clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
MvPSDEGs <- read.csv("Seurat/POA_IndependentAnalysis/POA_Fullmerge_30pcs_0.1_1.5cutoffMvP_SDEGS_Per_Clusterpct.csv", header=TRUE)
MvPSDEGs <- data.frame(MvPSDEGs)
MvUSDEGs <- read.csv("Seurat/POA_IndependentAnalysis/POA_Fullmerge_30pcs_0.1_1.5cutoffMvU_SDEGS_Per_Clusterpct.csv", header=TRUE)
MvUSDEGs <- data.frame(MvUSDEGs)
PvUEDGs <- read.csv("Seurat/POA_IndependentAnalysis/POA_Fullmerge_30pcs_0.1_1.5cutoffPvU_SDEGS_Per_Clusterpct.csv", header=TRUE)
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
FeaturePlot(mySeurat, features=c("MvP"), col=c("white","black"))
dev.off()
pdf(file=paste0(output,"_MvUSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("MvU"), col=c("white","black"))
dev.off()
pdf(file=paste0(output,"_PvUSDegsperclust.pdf"))
FeaturePlot(mySeurat, features=c("PvU"), col=c("white","black"))
dev.off() 
