library(Seurat)
library(ggplot2)
output <- "Seurat/VMH_IndependentAnalysis/VMH_Paper_UMAPheatMap_"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")

data <- read.table("VMHUMapTest.txt", header=TRUE, sep="\t")
head(data)

clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
clusters$MFrRatio <- data$ratio[match(clusters$clusters,data$id)]
head(clusters)
clustersMFr <- clusters["MFrRatio"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersMFr, col.name="MFrRatio")
head(mySeurat[[]])
pdf(file=paste0(output,"_VMHtaehongtest.pdf"))
FeaturePlot(mySeurat, features=c("MFrRatio"), cols=c("deeppink1", "blue"))
dev.off()

p1 <- FeaturePlot(mySeurat, features=c("MFrRatio"))
pdata <- p1$data
p <- ggplot(pdata, aes(x=UMAP_1, y=UMAP_2, color=MFrRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colours = c("deeppink1","gray","blue")) 
saveRDS(p1, file="TaehongTestUMAP.rds")
pdf(file=paste0(output,"VMHtaehongedited.pdf"))
p1
dev.off()