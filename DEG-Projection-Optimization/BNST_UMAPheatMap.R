library(Seurat)
library(ggplot2)
output <- "Seurat/BNST_IndependentAnalysis/BNST_Paper_UMAPheatMap_grayscale"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_20pcs.rds")

clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
MvPSDEGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_20pcs_0.1MvP_SDEGS_Per_Clusterpct.csv", header=TRUE)
MvPSDEGs <- data.frame(MvPSDEGs)
MvUSDEGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_20pcs_0.1MvU_SDEGS_Per_Clusterpct.csv", header=TRUE)
MvUSDEGs <- data.frame(MvUSDEGs)
PvUEDGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_20pcs_0.1PvU_SDEGS_Per_Clusterpct.csv", header=TRUE)
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
