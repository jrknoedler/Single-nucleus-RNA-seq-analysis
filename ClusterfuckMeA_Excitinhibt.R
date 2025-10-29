#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/MeA_IndependentAnalysis/MeA_Merged_Clusterfuck_Marchredo_0.1logcounts_"



all <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_CCAfiltered3_res1.2.rds")

mySeurat <- subset(all, idents=c("0","3","8","12","13","14","16","17","20","21","22","28","30","31","32","33"))

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/MeA_MvP_1.5cutoff.txt", header=FALSE)
eDEG <- read.table("Genelists/MeA_PvU_1.5cutoff.txt", header=FALSE)
sDEG2 <- read.table("Genelists/MeA_MvU_1.5cutoff.txt", header=FALSE)
eDEG <- unlist(eDEG)
sDEGs <- rbind(sDEG2, Dimorphic)
sDEGs <- unlist(sDEGs)
sDEGs <- unique(sDEGs)
sDEGs <- as.matrix(sDEGs)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered <- as.matrix(eDEG.filtered)

sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered <- as.matrix(sDEG2.filtered)
sDEGs.filtered <- intersect(genes.10x, sDEGs)
sDEGs.filtered <- as.matrix(sDEGs.filtered)
Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)
mySeurat <- NormalizeData(mySeurat)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)


avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=sDEGs.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in sDEGs.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "clusterspersdeg_excitatory.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "sdegspercluster_excitatory.csv"))

avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=eDEG.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in sDEGs.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "clustersperedeg_excitatory.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "edegspercluster_excitatory.csv"))






mySeurat <- subset(all, idents=c("0","3","8","12","13","14","16","17","20","21","22","28","30","31","32","33","25","4","9"), invert=TRUE)

genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/MeA_MvP_1.5cutoff.txt", header=FALSE)
eDEG <- read.table("Genelists/MeA_PvU_1.5cutoff.txt", header=FALSE)
sDEG2 <- read.table("Genelists/MeA_MvU_1.5cutoff.txt", header=FALSE)
eDEG <- unlist(eDEG)
sDEGs <- rbind(sDEG2, Dimorphic)
sDEGs <- unlist(sDEGs)
sDEGs <- unique(sDEGs)
sDEGs <- as.matrix(sDEGs)
Dimorphic <- as.matrix(Dimorphic)
eDEG <- as.matrix(eDEG)
sDEG2 <- as.matrix(sDEG2)
eDEG.filtered <- intersect(eDEG, genes.10x)
eDEG.filtered <- as.matrix(eDEG.filtered)

sDEG2.filtered <- intersect(sDEG2, genes.10x)
sDEG2.filtered <- as.matrix(sDEG2.filtered)
sDEGs.filtered <- intersect(genes.10x, sDEGs)
sDEGs.filtered <- as.matrix(sDEGs.filtered)
Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)
mySeurat <- NormalizeData(mySeurat)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)

avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=sDEGs.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in sDEGs.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "clusterspersdeg_inhibitory.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "sDEGs_Per_Cluster_inhibitory.csv"))

avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=eDEG.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in sDEGs.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.1 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "Clusters_Per_eDEG_inhibitory.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.1)
write.csv(counts, file=paste0(output, "eDEGs_Per_Cluster_inhibitory.csv"))
