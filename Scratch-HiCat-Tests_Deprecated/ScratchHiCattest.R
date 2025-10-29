#!/usr/bin/env Rscript



output <- "Scrattch_hicat/VMH/VMHMerged_filtered1_Paperanalysis_stricter"


library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)
library(Seurat)

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_filtered1_Paperanalysis.rds")
cells <- rownames(mySeurat[[]])
data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)
de.param <- de_param(padj.th = 0.05, lfc.th=0.5, low.th=1, q1.th=0.3, q2.th=NULL, q.diff.th=0.7, de.score.th=150, min.cells=10)
data <- cpm(data)
norm.data <- log2(data+1)
norm.data <- Matrix(norm.data, sparse=TRUE)
select.cells <- colnames(norm.data)
strict.param <- de_param(de.score.th=500)
onestep.result <- onestep_clust(norm.data, dim.method="pca", de.param=strict.param)
pdf(file=paste0(output,"_onestep.pdf"))
display.result <- display_cl(onestep.result$cl, norm.data, plot=TRUE, de.param=de.param)
display.result
dev.off()
iter.result <- iter_clust(norm.data, dim.method="pca", de.param=de.param, result=onestep.result)
pdf(file=paste0(output,"_iter.pdf"))
display_cl(iter.result$cl, norm.data, plot=TRUE, de.param=de.param)
dev.off()
rd.dat <- t(norm.data[iter.result$markers,])
merge.param <- de_param(de.score.th=200)
merge.result <- merge_cl(norm.data, cl=iter.result$cl, rd.dat=rd.dat, de.param=merge.param, return.markers=TRUE)
pdf(file=paste0(output,"_merge.pdf"))
display_cl(merge.result$cl, norm.data, plot=TRUE, de.param=merge.param)
dev.off()
iter.cl <- setNames(as.factor(iter.result$cl), select.cells)
iter.cl.df <- data.frame(cluster_id = unique(iter.cl), cluster_label=paste0("Pre-merge_cl_",unique(iter.cl)),cluster_color=rainbow(length(unique(iter.cl))))
rownames(iter.cl.df) <- iter.cl.df$cluster_id

pdf(file=paste0(output,"_compare.pdf"), height=20))
compare.result <- compare_annotate(merge.result$cl, iter.cl, iter.cl.df)
compare.result$g
dev.off()
cl <- compare.result$cl
cl.df <- compare.result$cl.df

pdf(file=paste0(output,"_DEgenes.pdf"))
display_cl(cl, norm.data, plot=TRUE, de.param=de.param, n.markers=10)
dev.off()
de.genes <- display.result$de.genes
#write.csv(de.genes, file=paste0(output,"_DEgenes.csv"))
set.seed(12345)
result <- run_consensus_clust(norm.data, niter=30, de.param=de.param, dim.method="pca", output_dir=paste0(output,"subsample_PCA"))
compare.result <- compare_annotate(result$cl.result$cl, iter.cl, iter.cl.df)
consensus.cl <- compare.result$cl
consensus.cl.df <- compare.result$cl.df
co.ratio <- result$co.result$co.ratio
#pdf(file=paste0(output,"_comatrixplot.pdf"))
#plot_co_matrix(co.ratio, consensus.cl, max.cl.size=50)
#dev.off()
cl.clean <- droplevels(consensus.cl)
head(cl.clean)
select.markers <- select_markers(norm.data, cl.clean, de.genes=de.genes, n.markers=50)
marker.genes <- select.markers$markers
de.genes <- select.markers$de.genes

cl.med <- ge_cl_medians(norm.data[marker.genes,], cl.clean)

1.rank <- setNames(1:nrow(consensus.cl.df), row.names(consensus.cl.df))

1.color <- setNames(as.character(consensus.cl.df$cluster_color), row.names(consensus.cl.df))

pdf(file=paste0(output,"_dendrogram.pdf"))
dend.result <- build_dend(cl.med[,levles(cl.clean)], 1.rank, 1.color, nboot=30)

dend <- dend.result$dendd

dend.labeled <- dend
labels(dend.labeled <- consensus.cl.df[lables(dend), "cluster_label"]
plot(dend.labeled)
dev.off()
pdf(file=paste0(output,"_dendroheatmap.pdf"))
cl.cor <- dend.results$cl.cor
row.names(cl.cor) <- colnames(cl.cor) <- consensus.cl.df[row.names(cl.cor), "cluster_label"]
heatmap.3(cl.cor, Rowv=dend, colv=dend, trace="non", col=heat.colors(100), cexRow=0.8, cexCol=0.8, breaks=c(-0.2, 0.2, seq(0.2, 1, length.out=99)))
dev.off()
