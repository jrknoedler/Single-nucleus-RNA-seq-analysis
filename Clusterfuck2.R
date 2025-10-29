#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(Matrix)
library(tidyverse)

output <- "Seurat/BNST_IndependentAnalysis/BNST_Merged_Clusterfuck_0.5logcounts_PrimedUPvMale"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_.rds")
genes.10x <- (x=rownames(x=mySeurat))
unlist(genes.10x)
Dimorphic <- read.table("Genelists/BNST_UpPrimedNirao.txt")
unlist(Dimorphic)
Dimorphic <- as.matrix(Dimorphic)
Dimorphic.filtered <- intersect(Dimorphic, genes.10x)
Dimorphic.filtered <- as.matrix(Dimorphic.filtered)
mySeurat <- NormalizeData(mySeurat)
data <- mySeurat@assays[["RNA"]]@data
data <- data.frame(data)
data <- t(data)
clusters <- mySeurat@active.ident
clusters <- data.frame(clusters)
merged <- cbind(data, clusters)
avg <- merged %>% group_by(clusters) %>% summarize_at(.vars=Dimorphic.filtered, .funs=c("mean"))
avg
avg %>% remove_rownames %>% column_to_rownames (var="clusters")
df_total= data.frame()
for (g in Dimorphic.filtered)
try({
results <- paste0(g,"results")
results = NULL
counts <- nrow(avg[avg[[g]] > 0.5 ,])
df <- data.frame(g, counts)
df_total = rbind(df_total, df)
})

write.csv(df_total, file=paste0(output, "Clusters_Per_SDEG.csv"))
avg <- data.matrix(avg)
counts <- rowSums(avg > 0.5)
write.csv(counts, file=paste0(output, "SDEGS_Per_Cluster.csv"))