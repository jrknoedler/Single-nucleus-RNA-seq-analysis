#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(scales)
output <- "Seurat/VMH_IndependentAnalysis/VMH_Merged_JuneFinal_10pct_bluescaletCTnorm"

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Fullmerge_Test_prepaperreclusterRound2_sexclude_malat1regress_filtered5_30pcsres1.2.rds")
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27")
names(new.cluster.ids) <- levels(mySeurat)
mySeurat <- RenameIdents(mySeurat, new.cluster.ids)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/VMH_MvP_1.5.txt")
eDEG <- read.table("Genelists/VMH_PvU_1.5.txt")
sDEG2 <- read.table("Genelists/VMH_MvU_1.5.txt")
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


AllsDegs <- read.table("topGO/Total_ByRegion/VMH_Genesonly.txt")
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
pct1 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 10), mvfrprop=MvFr/n())
head(pct1)
pct1 <- column_to_rownames(pct1, 'id')
head(pct1)



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


Pct3 <- DotPlot(mySeurat, features=Dimorphic.filtered)
data <- Pct3$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct3 <- data.df %>% group_by(id) %>% summarize(MfrTot=sum(pct.exp > 10))

pct3 <- column_to_rownames(pct3, 'id')





comb <- cbind(pct1, pct2, pct3)
head(comb)
MFrrat <- comb %>% summarize(MFrNorm=(MvFr-FrvM)/MfrTot)

MFrrat <- rownames_to_column(MFrrat, var="id")
head(MFrrat)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)

clusters$MFrNorm <- MFrrat$MFrNorm[match(clusters$clusters,MFrrat$id)]
clustersMFr <- clusters["MFrNorm"]
head(clustersMFr)
mySeurat <- AddMetaData(mySeurat, metadata=clustersMFr, col.name="MFrNorm")







Pct4 <- DotPlot(mySeurat, features=Up_MvFu.filtered)

data <- Pct4$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct4 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 10), mvfrprop=MvFr/n())
head(pct3)
pct4 <- column_to_rownames(pct4, 'id')
head(pct4)



Pct5 <- DotPlot(mySeurat, features=Up_FuvM.filtered)

data <- Pct5$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct5 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 10), frmprop=FrvM/n())

pct5 <- column_to_rownames(pct5, 'id')



Pct6 <- DotPlot(mySeurat, features=sDEG2.filtered)
data <- Pct6$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct6 <- data.df %>% group_by(id) %>% summarize(MfrTot=sum(pct.exp > 10))

pct6 <- column_to_rownames(pct6, 'id')



comb2 <- cbind(pct4, pct5, pct6)
comb2
MFurat <- comb2 %>% summarize(MFuNorm=(MvFr-FrvM)/MfrTot)
MFurat

MFurat <- rownames_to_column(MFurat, var="id")
head(MFurat)


clusters$MFuNorm <- MFurat$MFuNorm[match(clusters$clusters,MFurat$id)]
clustersMFu <- clusters["MFuNorm"]
head(clustersMFu)
mySeurat <- AddMetaData(mySeurat, metadata=clustersMFu, col.name="MFuNorm")





Pct7 <- DotPlot(mySeurat, features=Up_FrvFu.filtered)

data <- Pct7$data
head(data)
data
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct7 <- data.df %>% group_by(id) %>% summarize(MvFr=sum(pct.exp > 10), mvfrprop=MvFr/n())
head(pct5)
pct7 <- column_to_rownames(pct7, 'id')
head(pct7)



Pct8 <- DotPlot(mySeurat, features=Up_FuvFr.filtered)

data <- Pct8$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct8 <- data.df %>% group_by(id) %>% summarize(FrvM=sum(pct.exp > 10), frmprop=FrvM/n())
head(pct8)
pct8 <- column_to_rownames(pct8, 'id')
head(pct8)


Pct9 <- DotPlot(mySeurat, features=eDEG.filtered)
data <- Pct9$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct9 <- data.df %>% group_by(id) %>% summarize(MfrTot=sum(pct.exp > 10))

pct9 <- column_to_rownames(pct9, 'id')


comb <- cbind(pct7, pct8, pct9)
head(comb)
FrFurat <- comb %>% summarize(FrFuNorm=(MvFr-FrvM)/MfrTot)


FrFurat <- rownames_to_column(FrFurat, var="id")
head(FrFurat)


clusters$FrFuNorm <- FrFurat$FrFuNorm[match(clusters$clusters,FrFurat$id)]
clustersFrFu <- clusters["FrFuNorm"]
mySeurat <- AddMetaData(mySeurat, metadata=clustersFrFu, col.name="FrFuNorm")










p <- DimPlot(mySeurat, reduction="umap", label=TRUE, label.size=8)
pd <- p + theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(), axis.text=element_blank())+theme(legend.position="none")

blank <- ggplot()+theme_void()


p <- FeaturePlot(mySeurat, features=c("MFrNorm"))
pdata1 <- p$data
p1 <- ggplot(pdata1, aes(x=UMAP_1, y=UMAP_2, color=MFrNorm)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("deeppink1","gold","cyan","blue"),values= rescale(c(-1,0,0.5,1)),  guide="colorbar", limits=c(-1,1)) + theme(legend.title=element_blank(),legend.key.size=unit(1, "cm")) 
saveRDS(p1, file="TaehongTestUMAPauto.rds")
pdf(file=paste0(output,"VMHtaehongedited.pdf"))
p1
dev.off()

p <- FeaturePlot(mySeurat, features=c("MFuNorm"))
pdata2 <- p$data
p2 <- ggplot(pdata2, aes(x=UMAP_1, y=UMAP_2, color=MFuNorm)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("deeppink1","gold","cyan","blue"), values= rescale(c(-1,0,0.5,1)), guide="colorbar", limits=c(-1,1)) + theme(legend.title=element_blank(),  legend.key.size=unit(1, "cm"))

p <- FeaturePlot(mySeurat, features=c("FrFuNorm"))
pdata3 <- p$data
p3 <- ggplot(pdata3, aes(x=UMAP_1, y=UMAP_2, color=FrFuNorm)) + cowplot::theme_cowplot()+geom_point(size=0.3)+theme(axis.line=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), axis.title=element_blank()) + scale_colour_gradientn(colors=c("deeppink1","gold","cyan","blue"),values= rescale(c(-1,0,0.5,1)),  guide="colorbar", limits=c(-1,1)) + theme(legend.title=element_blank(),legend.key.size=unit(1, "cm"))

pdf(file=paste0(output,"Taehongcompiledbluescale.pdf"), width=30, height=10)
p1 + p2 + p3
dev.off()

pdf(file=paste0(output, "NewCombinedplotbluescale.pdf"), width=28)
(pd | blank| p1 |blank| p2 |blank| p3) + plot_layout(widths=c(1,0.5,1,0.5,1,0.5,1))
dev.off()
