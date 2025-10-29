#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_JuneFinal_25pct_"

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





p <- FeaturePlot(mySeurat, features=c("MFrRatio"))
pdata1 <- p$data
p1 <- ggplot(pdata1, aes(x=UMAP_1, y=UMAP_2, color=MFrRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colours = c("deeppink1","yellow","blue")) 



p <- FeaturePlot(mySeurat, features=c("MFuRatio"))
pdata2 <- p$data
p2 <- ggplot(pdata2, aes(x=UMAP_1, y=UMAP_2, color=MFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colours = c("dark green","yellow","blue")) 

p <- FeaturePlot(mySeurat, features=c("FrFuRatio"))
pdata3 <- p$data
p3 <- ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=FrFuRatio)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colours = c("dark green","yellow","deeppink1")) 

pdf(file=paste0(output,"Taehongcompiled.pdf"), width=30, height=10)
p1 + p2 + p3
dev.off()