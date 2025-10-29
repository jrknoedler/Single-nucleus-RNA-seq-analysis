#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(scales)
library(patchwork)
output <- "Seurat/BNST_IndependentAnalysis/BNST_UMAPtCTDEGlabels_JuneFinal_zratiotest"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)



Up_MvFr <- read.table("DEGmodules/BNST_UpMalevFr.txt", header=FALSE)
Up_MvFr <- unlist(Up_MvFr)
Up_MvFr.m <- as.matrix(Up_MvFr)
Up_MvFr.filtered <- intersect(Up_MvFr.m, genes.10x)

Up_FrvM <- read.table("DEGmodules/BNST_UpFrvMale.txt", header=FALSE)
Up_FrvM <- unlist(Up_FrvM)
Up_FrvM.m <- as.matrix(Up_FrvM)
Up_FrvM.filtered <- intersect(Up_FrvM.m, genes.10x)

Up_MvFu <- read.table("DEGmodules/BNST_UpMalevFu.txt", header=FALSE)
Up_MvFu <- unlist(Up_MvFu)
Up_MvFu.m <- as.matrix(Up_MvFu)
Up_MvFu.filtered <- intersect(Up_MvFu.m, genes.10x)

Up_FuvM <- read.table("DEGmodules/BNST_UpFuvMale.txt", header=FALSE)
Up_FuvM <- unlist(Up_FuvM)
Up_FuvM.m <- as.matrix(Up_FuvM)
Up_FuvM.filtered <- intersect(Up_FuvM.m, genes.10x)

Up_FrvFu <- read.table("DEGmodules/BNST_UpFrvFu.txt", header=FALSE)
Up_FrvFu <- unlist(Up_FrvFu)
Up_FrvFu.m <- as.matrix(Up_FrvFu)
Up_FrvFu.filtered <- intersect(Up_FrvFu.m, genes.10x)

Up_FuvFr <- read.table("DEGmodules/BNST_UpFuvFr.txt", header=FALSE)
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
pct1 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 10))
head(pct1)
pct1 <- column_to_rownames(pct1, 'id')
head(pct1)

pct1 <- scale(pct1)
pct1.df <- data.frame(pct1)
head(pct1)


Pct2 <- DotPlot(mySeurat, features=Up_FrvM.filtered)

data <- Pct2$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct2 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 10))
head(pct2)
pct2 <- column_to_rownames(pct2, 'id')
head(pct2)
pct2 <- scale(pct2)
pct2.df <- data.frame(pct2)
comb <- cbind(pct1.df, pct2.df)
head(comb)
MFrrat <- comb %>% summarize(MFrRatio=MvFr/FrvM)


MFrrat <- rownames_to_column(MFrrat, var="id")
head(MFrrat)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)

MZ <- rownames_to_column(pct1.df, var="id")
FZ <- rownames_to_column(pct2.df, var="id")


clusters$MFrRatio <- MFrrat$MFrRatio[match(clusters$clusters,MFrrat$id)]
clustersMFr <- clusters["MFrRatio"]
head(clustersMFr)
mySeurat <- AddMetaData(mySeurat, metadata=clustersMFr, col.name="MFrRatio")

clusters$MZS <- MZ$MvFr[match(clusters$clusters,MZ$id)]
clustersMZ<- clusters["MZS"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersMZ, col.name="MZ")
clusters$FZS <- FZ$FrvM[match(clusters$clusters,FZ$id)]
clustersFZ<- clusters["FZS"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersFZ, col.name="FZ")

p <- FeaturePlot(mySeurat, features=c("MFrRatio"))
pdata1 <- p$data
p1 <- ggplot(pdata1, aes(x=UMAP_1, y=UMAP_2, color=MFrRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("deeppink1","gold","cyan","blue"), guide="colorbar") + theme(legend.title=element_blank(), legend.key.size=unit(1, "cm"))
saveRDS(p1, file="TaehongTestUMAPauto.rds")
pdf(file=paste0(output,"VMHtaehongedited.pdf"))
p1
dev.off()

p <- FeaturePlot(mySeurat, features=c("MZ"))
pdata2 <- p$data
p2 <- ggplot(pdata2, aes(x=UMAP_1, y=UMAP_2, color=MZ)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("deeppink1","gold","cyan","blue"),  guide="colorbar") + theme(legend.title=element_blank(), legend.key.size=unit(1, "cm"))

p <- FeaturePlot(mySeurat, features=c("FZ"))
pdata3 <- p$data
p3 <- ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=FZ)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("deeppink1","gold","cyan","blue"),  guide="colorbar") + theme(legend.title=element_blank(), legend.key.size=unit(1, "cm"))

pdf(file=paste0(output,"zscorestest.pdf"), width=30, height=10)
p1 + p2 + p3
dev.off()