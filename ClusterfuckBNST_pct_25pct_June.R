#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)


output <- "Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_Juneredo_25pct_"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound4_sexclude_malat1regress_filtered5_2res1.2_30pcsredo.rds")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/BNST_MvP_1.5cutoff.txt")
eDEG <- read.table("Genelists/BNST_PvU_1.5cutoff.txt")
sDEG2 <- read.table("Genelists/BNST_MvU_1.5cutoff.txt")
unlist(eDEG)
unlist(sDEG2)
unlist(Dimorphic)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered <- as.matrix(eDEG.filtered)

sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered <- as.matrix(sDEG2.filtered)

Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)


AllsDegs <- read.table("topGO/Total_ByRegion/BNST_Genesonly.txt")
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


mySeurat <- NormalizeData(mySeurat)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)


Pct <- DotPlot(mySeurat, features=Dimorphic.filtered)
data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 25), mvpprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"MvP_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), mvpprop=num/n_clust)
write.csv(dist, file=paste0(output,"MvP_clustersperDEG.csv"))
dist



Pct <- DotPlot(mySeurat, features=eDEG.filtered)
data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 25), pvuprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"PvU_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), pvuprop=num/n_clust)
head(dist)
write.csv(dist, file=paste0(output,"PvU_clustersperDEG.csv"))

Pct <- DotPlot(mySeurat, features=sDEG2.filtered)
data <- Pct$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
pct <- data.df %>% group_by(id) %>% summarize(n_cells=n(), num=sum(pct.exp > 25), mvuprop=num/n_cells)
pct
write.csv(pct, file=paste0(output,"MvU_sDEGspercluster.csv"))
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >25), mvuprop=num/n_clust)
head(dist)
write.csv(dist, file=paste0(output,"MvU_clustersperDEG.csv"))
allclust <- dist[dist$num == 35,]
allclust

genestotal=data.frame()
for (g in AllsDegs.filtered)
try({

pct <- sum(GetAssayData(mySeurat, slot="data")[g,]>0)/nrow(mySeurat@meta.data)
df <- data.frame(g, pct)
genestotal=rbind(genestotal, df)
})
write.csv(genestotal, file=paste0(output,"allDEGs_allcellspct.csv"))
Bigplot <- DotPlot(mySeurat, features=AllsDegs.filtered)

data <- Bigplot$data
head(data)
data.df <- data.frame(data)
data.df <- rownames_to_column(data.df, var="gene")
head(data.df)
dist <- data.df %>% group_by(features.plot) %>% summarize(n_clust=n(), num=sum(pct.exp >5), totalprop=num/n_clust)
head(dist)
allclust <- dist[dist$num == 35,]
allclust
write.csv(allclust, file=paste0(output,"allubiquitousDEGs.csv"))
unique <- dist[dist$num == 1, ]
unique
uniqueclust <- unique$features.plot
uniqueclust.df <- data.frame(uniqueclust)
uniqueclustgenes <- data.df %>% filter(features.plot %in% uniqueclust)
uniqueclustassign <- uniqueclustgenes[uniqueclustgenes$pct.exp > 5,]
uniqueclustassign
write.csv(uniqueclustassign, file=paste0(output,"BNST_clusterspecDEGs.csv"))