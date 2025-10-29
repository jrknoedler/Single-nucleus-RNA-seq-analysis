#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(scales)
library(patchwork)
output <- "X:/Seurat/BNST_IndependentAnalysis/BNST_UMAPtCTDEGlabels_JuneFinal_10pct_zscaleratio"

mySeurat <- readRDS("X:/Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("X:/Genelists/BNST_MvP_1.5cutoff.txt")
eDEG <- read.table("X:/Genelists/BNST_PvU_1.5cutoff.txt")
sDEG2 <- read.table("X:/Genelists/BNST_MvU_1.5cutoff.txt")
unlist(eDEG)
unlist(sDEG2)
unlist(Dimorphic)
allSDEGs <- union(Dimorphic, sDEG2)
allSDEGs <- as.matrix(allSDEGs)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered <- as.matrix(eDEG.filtered)
allSDEGs.filtered <- intersect(allSDEGs, genes.10x)
allSDEGs.filtered <- as.matrix(allSDEGs.filtered)
sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered <- as.matrix(sDEG2.filtered)

Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)


AllsDegs <- read.table("X:/topGO/Total_ByRegion/BNST_Genesonly.txt")
AllsDegs <- unlist(AllsDegs)
AllsDegs <- as.matrix(AllsDegs)
AllsDegs.filtered <- intersect(AllsDegs, genes.10x)
AllsDegs.filtered <- as.matrix(AllsDegs.filtered)
totdegs <- nrow(AllsDegs)
detecdegs <- nrow(AllsDegs.filtered)
conc <- cbind(detecdegs, totdegs)
conc.df <- data.frame(conc)
write.csv(conc.df, file=paste0(output,"10xtrapconcordance.csv"))
AllsDegs.filtered <- as.matrix(AllsDegs.filtered)



Up_MvFr <- read.table("X:/DEGmodules/BNST_UpMalevFr.txt", header=FALSE)
Up_MvFr <- unlist(Up_MvFr)
Up_MvFr.m <- as.matrix(Up_MvFr)
Up_MvFr.filtered <- intersect(Up_MvFr.m, genes.10x)

Up_FrvM <- read.table("X:/DEGmodules/BNST_UpFrvMale.txt", header=FALSE)
Up_FrvM <- unlist(Up_FrvM)
Up_FrvM.m <- as.matrix(Up_FrvM)
Up_FrvM.filtered <- intersect(Up_FrvM.m, genes.10x)

Up_MvFu <- read.table("X:/DEGmodules/BNST_UpMalevFu.txt", header=FALSE)
Up_MvFu <- unlist(Up_MvFu)
Up_MvFu.m <- as.matrix(Up_MvFu)
Up_MvFu.filtered <- intersect(Up_MvFu.m, genes.10x)

Up_FuvM <- read.table("X:/DEGmodules/BNST_UpFuvMale.txt", header=FALSE)
Up_FuvM <- unlist(Up_FuvM)
Up_FuvM.m <- as.matrix(Up_FuvM)
Up_FuvM.filtered <- intersect(Up_FuvM.m, genes.10x)

Up_FrvFu <- read.table("X:/DEGmodules/BNST_UpFrvFu.txt", header=FALSE)
Up_FrvFu <- unlist(Up_FrvFu)
Up_FrvFu.m <- as.matrix(Up_FrvFu)
Up_FrvFu.filtered <- intersect(Up_FrvFu.m, genes.10x)

Up_FuvFr <- read.table("X:/DEGmodules/BNST_UpFuvFr.txt", header=FALSE)
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



Pct1 <- DotPlot(mySeurat, features=Up_MvFr.filtered)

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



Pct2 <- DotPlot(mySeurat, features=Up_FrvM.filtered)

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
mySeurat <- AddMetaData(mySeurat, metadata=clustersMFr, col.name="MFrRatio")



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

pdf(file=paste0(output,"_VMHtaehongtestauto.pdf"))
FeaturePlot(mySeurat, features=c("MFrRatio"), cols=c("deeppink1", "blue"))
dev.off()

p <- DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=5)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")

blank <- ggplot()+theme_void()

p <- FeaturePlot(mySeurat, features=c("AllDEGs"))
pgdat <- p$data
pg <- ggplot(pgdat, aes(x=UMAP_1, y=UMAP_2, color=AllDEGs)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("white","black"),   guide="colorbar", limits=c(0.2,0.8)) + theme(legend.position="none") 

blank <- ggplot()+theme_void()

pt <- DimPlot(mySeurat, reduction="umap", label=FALSE, cols=c("1"="orangered1","2"="orangered1","3"="orangered1","4"="orangered1","5"="orangered1","6"="light blue","7"="orangered1","8"="orangered1","9"="orangered1","10"="orangered1","11"="orangered1","12"="orangered1","13"="orangered1","14"="orangered1","15"="light blue","16"="orangered1","17"="orangered1","18"="orangered1","19"="orangered1","20"="light blue","21"="orangered1","22"="orangered1","23"="orangered1","24"="orangered1","25"="orangered1","26"="orangered1","27"="orangered1","28"="orangered1","29"="orangered1","30"="orangered1","31"="orangered1","32"="orangered1","33"="orangered1","34"="orangered1","35"="light blue","36"="orangered1"))
ptf <- pt + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")


p <- FeaturePlot(mySeurat, features=c("MFrRatio"))
pdata1 <- p$data
p1 <- ggplot(pdata1, aes(x=UMAP_1, y=UMAP_2, color=MFrRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradient2(low="deeppink1", high="cyan3", mid="light gray", na.value="royalblue3", guide="colorbar", limits=c(-2,2)) + theme(legend.position="none") 
saveRDS(p1, file="TaehongTestUMAPauto.rds")
pdf(file=paste0(output,"VMHtaehongedited.pdf"))
p1
dev.off()

p <- FeaturePlot(mySeurat, features=c("MFuRatio"))
pdata2 <- p$data
p2 <- ggplot(pdata2, aes(x=UMAP_1, y=UMAP_2, color=MFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradient2(low="deeppink1", high="cyan3", mid="light gray", na.value="royalblue3", guide="colorbar",limits=c(-2,2)) + theme(legend.position="none") 

p <- FeaturePlot(mySeurat, features=c("FrFuRatio"))
pdata3 <- p$data
p3 <- ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=FrFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradient2(low="deeppink1", high="cyan3", mid="light gray", na.value="royalblue3", guide="colorbar", limits=c(-2,2))+ theme(legend.position="none") 

pdf(file=paste0(output,"Taehongcompiledbluescale.pdf"), width=10, height=40)
p1 + p2 + p3
dev.off()
pc <- DimPlot(mySeurat, reduction="umap", label=FALSE, cols=c("royalblue3", "deeppink1","chartreuse3"), group.by="Hormone", shuffle=TRUE)
pce <- pc + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank(), plot.title=element_blank())+theme(legend.position="none")

pdf(file=paste0(output,"celltyperep_R5.pdf"))
pce
dev.off()
pdf(file=paste0(output, "NewCombinedplotbluescalevertsixstack.pdf"), height=24, width=3.5)
(pd/pce/ptf/pg/p1/p2/p3) 
dev.off()