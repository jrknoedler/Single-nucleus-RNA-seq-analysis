library(Seurat)
library(ggplot2)
library(patchwork)
output <- "Seurat/BNST_IndependentAnalysis/BNST_Paper_UMAPheatMap_25pct_grayscale_final"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")

clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
MvPSDEGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Juneredo_25pct_MvP_sDEGspercluster.csv", header=TRUE)
MvPSDEGs <- data.frame(MvPSDEGs)
head(MvPSDEGs)
MvUSDEGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Juneredo_25pct_MvU_sDEGspercluster.csv", header=TRUE)
MvUSDEGs <- data.frame(MvUSDEGs)
PvUEDGs <- read.csv("Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Juneredo_25pct_PvU_sDEGspercluster.csv", header=TRUE)
PvUEDGs <- data.frame(PvUEDGs)
clusters$MvP <- MvPSDEGs$mvpprop[match(clusters$clusters,MvPSDEGs$id)]
head(clusters)
clustersMvP <- clusters["MvP"]
clusters$MvU <- MvUSDEGs$mvuprop[match(clusters$clusters,MvUSDEGs$id)]
clustersMvU <- clusters["MvU"]
clusters$PvU <- PvUEDGs$pvuprop[match(clusters$clusters,PvUEDGs$id)]
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
FeaturePlot(mySeurat, features=c("PvU"), col=c("white","black"))
dev.off()

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)

pMvP <- FeaturePlot(mySeurat, features=c("MvP"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
pMvU <- FeaturePlot(mySeurat, features=c("MvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)
pPvU <- FeaturePlot(mySeurat, features=c("PvU"), col=c("white","black"), min.cutoff=0, max.cutoff=1)

pMvPdat <- pMvP$data
pMvUdat <- pMvU$data
pPvUdat <- pPvU$data

p <- DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=8)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")
p1 <- ggplot(pMvPdat, aes(x=UMAP_1, y=UMAP_2, color=MvP)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) +theme(legend.position="none") 



p2 <- ggplot(pMvUdat, aes(x=UMAP_1, y=UMAP_2, color=MvU))+ cowplot::theme_cowplot() +geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) +theme(legend.position="none") 

blank <- ggplot()+theme_void()

p3 <- ggplot(pPvUdat, aes(x=UMAP_1, y=UMAP_2, color=PvU))+ cowplot::theme_cowplot() +geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colours = c("white","black"), limits=c(0.1, 0.7)) +theme(legend.position="none") 
pdf(file=paste0(output, "Combinedplot.pdf"), width=28)
(pd | blank| p1 |blank| p2 |blank| p3) + plot_layout(widths=c(1,0.5,1,0.5,1,0.5,1))
dev.off()