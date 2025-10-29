#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(scales)
memory.limit(size=64000)
output <- "X:/Seurat/MeA_IndependentAnalysis/MeA_Merged_DEGratios_JuneFinal_10pct_ratiozscaled_vert"

mySeurat <- readRDS("X:/Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered1_res1.2marredo_esr1filt_33filt2_1.2_30pcs.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)

mySeuratx <- subset(mySeurat, idents=c("27"), invert=TRUE)
mySeuraty <- subset(mySeurat, idents=c("27"))
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)


AllsDegs <- read.table("X:/topGO/Total_ByRegion/MeA_Genesonly.txt")
AllsDegs <- unlist(AllsDegs)
AllsDegs <- as.matrix(AllsDegs)
AllsDegs.filtered <- intersect(AllsDegs, genes.10x)
AllsDegs.filtered <- as.matrix(AllsDegs.filtered)

Up_MvFr <- read.table("X:/DEGmodules/MeA_UpMvFr.txt", header=FALSE)
Up_MvFr <- unlist(Up_MvFr)
Up_MvFr.m <- as.matrix(Up_MvFr)
Up_MvFr.filtered <- intersect(Up_MvFr.m, genes.10x)

Up_FrvM <- read.table("X:/DEGmodules/MeA_UpFrM.txt", header=FALSE)
Up_FrvM <- unlist(Up_FrvM)
Up_FrvM.m <- as.matrix(Up_FrvM)
Up_FrvM.filtered <- intersect(Up_FrvM.m, genes.10x)

Up_MvFu <- read.table("X:/DEGmodules/MeA_UpMvFu.txt", header=FALSE)
Up_MvFu <- unlist(Up_MvFu)
Up_MvFu.m <- as.matrix(Up_MvFu)
Up_MvFu.filtered <- intersect(Up_MvFu.m, genes.10x)

Up_FuvM <- read.table("X:/DEGmodules/MeA_UpFuvM.txt", header=FALSE)
Up_FuvM <- unlist(Up_FuvM)
Up_FuvM.m <- as.matrix(Up_FuvM)
Up_FuvM.filtered <- intersect(Up_FuvM.m, genes.10x)

Up_FrvFu <- read.table("X:/DEGmodules/MeA_UpFrvFu.txt", header=FALSE)
Up_FrvFu <- unlist(Up_FrvFu)
Up_FrvFu.m <- as.matrix(Up_FrvFu)
Up_FrvFu.filtered <- intersect(Up_FrvFu.m, genes.10x)

Up_FuvFr <- read.table("X:/DEGmodules/MeA_UpFuvFr.txt", header=FALSE)
Up_FuvFr <- unlist(Up_FuvFr)
Up_FuvFr.m <- as.matrix(Up_FuvFr)
Up_FuvFr.filtered <- intersect(Up_FuvFr.m, genes.10x)

data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)






Pcta <- DotPlot(mySeurat, features=AllsDegs.filtered)
data <- Pcta$data

data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pcta <- data.df %>% group_by(id) %>% summarize(AllDEGs=sum(pct.exp > 10), DEGprop=AllDEGs/n())

clusters$DEGprop <- pcta$DEGprop[match(clusters$clusters,pcta$id)]
clustersall <- clusters["DEGprop"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersall, col.name="AllDEGs")





Pct1 <- DotPlot(mySeuratx, features=Up_MvFr.filtered)

data <- Pct1$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct1 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 10), mvfrprop=MvFr/n())
head(pct1)
pct1 <- column_to_rownames(pct1, 'id')
head(pct1)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_MvFrclustersperDEG.csv"))




Pct2 <- DotPlot(mySeuratx, features=Up_FrvM.filtered)

data <- Pct2$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct2 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 10), frmprop=FrvM/n())
head(pct2)
pct2 <- column_to_rownames(pct2, 'id')
head(pct2)
write.csv(pct2, file=paste0(output,"FroverM_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_FrvMclustersperDEG.csv"))



comb <- cbind(pct1, pct2)
head(comb)
MFrrat <- comb %>% summarize(MFrRatio=mvfrprop/frmprop)


MFrrat <- scale(MFrrat)
MFrrat <- data.frame(MFrrat)

MFrrat <- rownames_to_column(MFrrat, var="id")
head(MFrrat)


clusters$MFrRatio <- MFrrat$MFrRatio[match(clusters$clusters,MFrrat$id)]
clustersMFr <- clusters["MFrRatio"]
head(clustersMFr)
mySeuratx <- AddMetaData(mySeuratx, metadata=clustersMFr, col.name="MFrRatio")



Pct3 <- DotPlot(mySeurat, features=Up_MvFu.filtered)

data <- Pct3$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct3 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 10), mvfrprop=MvFr/n())
head(pct3)
pct3 <- column_to_rownames(pct3, 'id')
head(pct3)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_MvFuclustersperDEG.csv"))




Pct4 <- DotPlot(mySeurat, features=Up_FuvM.filtered)

data <- Pct4$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct4 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 10), frmprop=FrvM/n())
head(pct4)
pct4 <- column_to_rownames(pct4, 'id')
head(pct4)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_FuvMclustersperDEG.csv"))



comb2 <- cbind(pct3, pct4)
comb2
MFurat <- comb2 %>% summarize(MFuRatio=mvfrprop/frmprop)
MFurat

MFurat <- scale(MFurat)
MFurat <- data.frame(MFurat)


MFurat <- rownames_to_column(MFurat, var="id")
head(MFurat)


clusters$MFuRatio <- MFurat$MFuRatio[match(clusters$clusters,MFurat$id)]
clustersMFu <- clusters["MFuRatio"]
head(clustersMFu)
mySeurat <- AddMetaData(mySeurat, metadata=clustersMFu, col.name="MFuRatio")



Pct5 <- DotPlot(mySeurat, features=Up_FrvFu.filtered)

data <- Pct5$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct5 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 10), mvfrprop=MvFr/n())
head(pct5)
pct5 <- column_to_rownames(pct5, 'id')
head(pct5)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_FrvFuclustersperDEG.csv"))




Pct6 <- DotPlot(mySeurat, features=Up_FuvFr.filtered)

data <- Pct6$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct6 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 10), frmprop=FrvM/n())
head(pct6)
pct6 <- column_to_rownames(pct6, 'id')
head(pct6)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >10), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"Up_FuvFrclustersperDEG.csv"))

p2t <- DimPlot(mySeurat, reduction="umap", label=FALSE, cols=c("1"="light blue","2"="orangered1","3"="orangered1","4"="light blue","5"="orangered1","6"="orangered1","7"="light blue","8"="orangered1","9"="light blue","10"="light blue","11"="orangered1","12"="orangered1","13"="light blue","14"="light blue","15"="light blue","16"="light blue","17"="orangered1","18"="light blue","19"="orangered1","20"="orangered1","21"="light blue","22"="light blue","23"="light blue","24"="light blue","25"="orangered1","26"="orangered1","27"="orangered1","28"="orangered1","29"="light blue","30"="light blue","31"="light blue","32"="light blue","33"="light blue","34"="orangered1"))
p2f <- p2t +theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")


comb <- cbind(pct5, pct6)
head(comb)
FrFurat <- comb %>% summarize(FrFuRatio=mvfrprop/frmprop)

FrFurat <- scale(FrFurat)
FrFurat <- data.frame(FrFurat)



FrFurat <- rownames_to_column(FrFurat, var="id")
head(FrFurat)


clusters$FrFuRatio <- FrFurat$FrFuRatio[match(clusters$clusters,FrFurat$id)]
clustersFrFu <- clusters["FrFuRatio"]
head(clustersMFu)
mySeurat <- AddMetaData(mySeurat, metadata=clustersFrFu, col.name="FrFuRatio")


breaks1 <- c(0.3,1,2)

p <- DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=5)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")

blank <- ggplot()+theme_void()

p <- FeaturePlot(mySeurat, features=c("AllDEGs"))
pgdat <- p$data
pg <- ggplot(pgdat, aes(x=UMAP_1, y=UMAP_2, color=AllDEGs)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("white","black"),   guide="colorbar", limits=c(0.2,0.8)) + theme(legend.position="none")

blank <- ggplot()+theme_void()


head(mySeuratx[[]])
p <- FeaturePlot(mySeuratx, features=c("MFrRatio"))
pdata1 <- p$data
pprime <- DimPlot(mySeuraty, reduction="umap")
pprimedat <- pprime$data
p1 <- ggplot(pdata1, aes(x=UMAP_1, y=UMAP_2, color=MFrRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradient2(low="deeppink1", high="cyan3", mid="light gray", na.value="blue", guide="colorbar", limits=c(-2,2), breaks=c(0))  + theme(legend.position="none") + geom_point(data=pprimedat, aes(x=UMAP_1,y=UMAP_2), size=0.3, color="light gray")+cowplot::theme_cowplot() + theme(axis.line=element_blank(), axis.text=element_blank()) +theme(axis.title=element_blank(), axis.ticks=element_blank())  +theme(legend.position="none")
saveRDS(p1, file="TaehongTestUMAPauto.rds")
pdf(file=paste0(output,"VMHtaehongedited.pdf"))
p1
dev.off()

p <- FeaturePlot(mySeurat, features=c("MFuRatio"))
pdata2 <- p$data
p2 <- ggplot(pdata2, aes(x=UMAP_1, y=UMAP_2, color=MFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradient2(low="deeppink1", high="cyan3", mid="light gray", na.value="blue", guide="colorbar", limits=c(-2,2), breaks=c(0))  + theme(legend.position="none")

p <- FeaturePlot(mySeurat, features=c("FrFuRatio"))
pdata3 <- p$data
p3 <- ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=FrFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradient2(low="deeppink1", high="cyan3", mid="light gray", na.value="blue", guide="colorbar", limits=c(-2,2), breaks=c(0))  + theme(legend.position="none")

pdf(file=paste0(output,"Taehongcompiledbluescale.pdf"), width=10, height=40)
p1 + p2 + p3
dev.off()

pdf(file=paste0(output, "NewCombinedplotbluescalevert.pdf"), height=21, width=3.5)
(pd /p2f/pg/  p1 / p2 / p3) 
dev.off()