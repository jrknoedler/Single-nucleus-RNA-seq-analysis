#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(scales)
output <- "Seurat/VMH_IndependentAnalysis/VMH_Merged_JuneFinal_25pct_bluescale"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)



Up_MvFr <- read.table("DEGmodules/VMH_MaleupvPrimed.txt", header=FALSE)
Up_MvFr <- unlist(Up_MvFr)
Up_MvFr.m <- as.matrix(Up_MvFr)
Up_MvFr.filtered <- intersect(Up_MvFr.m, genes.10x)

Up_FrvM <- read.table("DEGmodules/VMH_UpPrimedvMale.txt", header=FALSE)
Up_FrvM <- unlist(Up_FrvM)
Up_FrvM.m <- as.matrix(Up_FrvM)
Up_FrvM.filtered <- intersect(Up_FrvM.m, genes.10x)

Up_MvFu <- read.table("DEGmodules/VMH_UpMalevUnprimed.txt", header=FALSE)
Up_MvFu <- unlist(Up_MvFu)
Up_MvFu.m <- as.matrix(Up_MvFu)
Up_MvFu.filtered <- intersect(Up_MvFu.m, genes.10x)

Up_FuvM <- read.table("DEGmodules/VMH_UpUnprimedvMale.txt", header=FALSE)
Up_FuvM <- unlist(Up_FuvM)
Up_FuvM.m <- as.matrix(Up_FuvM)
Up_FuvM.filtered <- intersect(Up_FuvM.m, genes.10x)

Up_FrvFu <- read.table("DEGmodules/VMH_UpPrimedvFu.txt", header=FALSE)
Up_FrvFu <- unlist(Up_FrvFu)
Up_FrvFu.m <- as.matrix(Up_FrvFu)
Up_FrvFu.filtered <- intersect(Up_FrvFu.m, genes.10x)

Up_FuvFr <- read.table("DEGmodules/VMH_UpFuvPrimed.txt", header=FALSE)
Up_FuvFr <- unlist(Up_FuvFr)
Up_FuvFr.m <- as.matrix(Up_FuvFr)
Up_FuvFr.filtered <- intersect(Up_FuvFr.m, genes.10x)

data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)


Pct1 <- DotPlot(mySeurat, features=Up_MvFr.filtered)

data <- Pct1$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct1 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 25), mvfrprop=MvFr/n())
head(pct1)
pct1 <- column_to_rownames(pct1, 'id')
head(pct1)



Pct2 <- DotPlot(mySeurat, features=Up_FrvM.filtered)

data <- Pct2$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct2 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 25), frmprop=FrvM/n())
head(pct2)
pct2 <- column_to_rownames(pct2, 'id')
head(pct2)
write.csv(pct2, file=paste0(output,"FroverM_sDEGspercluster.csv"))

comb <- cbind(pct1, pct2)
head(comb)
MFrrat <- comb %>% summarize(MFrRatio=mvfrprop/frmprop)

MFrrat <- rownames_to_column(MFrrat, var="id")
head(MFrrat)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)

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
pct3 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 25), mvfrprop=MvFr/n())
head(pct3)
pct3 <- column_to_rownames(pct3, 'id')
head(pct3)



Pct4 <- DotPlot(mySeurat, features=Up_FuvM.filtered)

data <- Pct4$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct4 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 25), frmprop=FrvM/n())
head(pct4)
pct4 <- column_to_rownames(pct4, 'id')
head(pct4)

comb2 <- cbind(pct3, pct4)
comb2
MFurat <- comb2 %>% summarize(MFuRatio=mvfrprop/frmprop)
MFurat

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
pct5 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 25), mvfrprop=MvFr/n())
head(pct5)
pct5 <- column_to_rownames(pct5, 'id')
head(pct5)



Pct6 <- DotPlot(mySeurat, features=Up_FuvFr.filtered)

data <- Pct6$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct6 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 25), frmprop=FrvM/n())
head(pct6)
pct6 <- column_to_rownames(pct6, 'id')
head(pct6)

comb <- cbind(pct5, pct6)
head(comb)
FrFurat <- comb %>% summarize(FrFuRatio=mvfrprop/frmprop)


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

p <- DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=8)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")

blank <- ggplot()+theme_void()


p <- FeaturePlot(mySeurat, features=c("MFrRatio"))
pdata1 <- p$data
p1 <- ggplot(pdata1, aes(x=UMAP_1, y=UMAP_2, color=MFrRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("dark green","gold","deeppink1"), values=rescale(c(0,1,1.8)),  guide="colorbar",limits=c(0,1.8), breaks=c(0,1,1.8)) + theme(legend.title=element_blank(), legend.text=element_blank(), legend.key.size=unit(1, "cm")) 
saveRDS(p1, file="TaehongTestUMAPauto.rds")
pdf(file=paste0(output,"VMHtaehongedited.pdf"))
p1
dev.off()

p <- FeaturePlot(mySeurat, features=c("MFuRatio"))
pdata2 <- p$data
p2 <- ggplot(pdata2, aes(x=UMAP_1, y=UMAP_2, color=MFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("dark green","gold","deeppink1"),  values=rescale(c(0,1,1.7)), guide="colorbar", limits=c(0,1.7), breaks=c(0,1,1.7)) + theme(legend.title=element_blank(), legend.text=element_blank(), legend.key.size=unit(1, "cm"))

p <- FeaturePlot(mySeurat, features=c("FrFuRatio"))
pdata3 <- p$data
p3 <- ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=FrFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("dark green","gold","deeppink1"),values=rescale(c(0,1,3)),  guide="colorbar",limits=c(0,3), breaks=c(0,1,3)) + theme(legend.title=element_blank(), legend.text=element_blank(), legend.key.size=unit(1, "cm"))

pdf(file=paste0(output,"Taehongcompiledredscaled.pdf"), width=30, height=10)
p1 + p2 + p3
dev.off()

pdf(file=paste0(output, "NewCombinedplotredscaled.pdf"), width=28)
(pd | blank| p1 |blank| p2 |blank| p3) + plot_layout(widths=c(1,0.5,1,0.5,1,0.5,1))
dev.off()