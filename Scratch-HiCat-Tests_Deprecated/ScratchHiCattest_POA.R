#!/usr/bin/env Rscript



output <- "Scrattch_hicat/POA/POAMerged_Paperanalysis_Sun"


library(scrattch.hicat)
library(dendextend)
library(dplyr)
library(matrixStats)
library(Matrix)
library(Seurat)

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Test.rds")
cells <- rownames(mySeurat[[]])
data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)
de.param <- de_param(padj.th = 0.05, lfc.th=0.6, low.th=1, q1.th=0.3, q2.th=NULL, q.diff.th=0.1, de.score.th=100, min.cells=10)
data <- cpm(data)
norm.data <- log2(data+1)
norm.data <- Matrix(norm.data, sparse=TRUE)
select.cells <- colnames(norm.data)
gene.counts <- colSums(norm.data > 0)
rm.eigen <- matrix(log2(gene.counts), ncol=1)
row.names(rm.eigen) <- names(gene.counts)
colnames(rm.eigen) <- "log2GeneCounts"
strict.param <- de_param(de.score.th=500)
strict.param <- de_param(q.diff.th=0.9)
onestep.result <- onestep_clust(norm.data, dim.method="pca", max.dim=30, de.param=strict.param, rm.eigen=rm.eigen)
pdf(file=paste0(output,"_onestep.pdf"))
display.result <- display_cl(onestep.result$cl, norm.data, plot=TRUE, de.param=strict.param, n.markers=5)
display.result
dev.off()
saveRDS(onestep.result, file=paste0(output,"_strictparam.rds"))
iter.result <- iter_clust(norm.data, dim.method="pca", de.param=de.param, result=onestep.result, rm.eigen=rm.eigen)
saveRDS(iter.result, file=paste0(output,"_interatedclust.rds"))
iter.markers <- iter.result$markers
write.csv(iter.markers, file=paste0(output,"_interatedmarkers.csv"))
pdf(file=paste0(output,"_iter.pdf"))
display_cl(iter.result$cl, norm.data, plot=TRUE, de.param=de.param, n.markers=5)
dev.off()
rd.dat <- t(norm.data[iter.result$markers,])
merge.param1 <- de_param(de.score.th=150)
merge.result1 <- merge_cl(norm.data, cl=iter.result$cl, rd.dat=rd.dat, de.param=merge.param1, return.markers=TRUE)
saveRDS(merge.result1, file=paste0(output,"_Merged1.rds"))
pdf(file=paste0(output,"_merge1.pdf"))
display_cl(merge.result1$cl, norm.data, plot=TRUE, de.param=merge.param1, n.markers=5)
dev.off()
merge.param2 <- de_param(de.score.th=300)
merge.result2 <- merge_cl(norm.data, cl=iter.result$cl, rd.dat=rd.dat, de.param=merge.param2, return.markers=TRUE)
saveRDS(merge.result2, file=paste0(output,"_Merged2.rds"))
pdf(file=paste0(output,"_merge2.pdf"))
display_cl(merge.result2$cl, norm.data, plot=TRUE, de.param=merge.param2, n.markers=5)
dev.off()
iter.cl <- setNames(as.factor(iter.result$cl), select.cells)
iter.cl.df <- data.frame(cluster_id = unique(iter.cl), cluster_label=paste0("Pre-merge_cl_",unique(iter.cl)),cluster_color=rainbow(length(unique(iter.cl))))
rownames(iter.cl.df) <- iter.cl.df$cluster_id

pdf(file=paste0(output,"_compare1.pdf"), height=20)
compare.result1 <- compare_annotate(merge.result1$cl, iter.cl, iter.cl.df)
compare.result1$g
dev.off()
pdf(file=paste0(output,"_compare2.pdf"), height=20)
compare.result2 <- compare_annotate(merge.result2$cl, iter.cl, iter.cl.df)
compare.result2$g
dev.off()
cl <- compare.result1$cl
cl.df <- compare.result1$cl.df

pdf(file=paste0(output,"_DEgenes.pdf"))
display_cl(cl, norm.data, plot=TRUE, de.param=de.param, n.markers=5)
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

cl.med <- get_cl_medians(norm.data[marker.genes,], cl.clean)

rank <- setNames(1:nrow(consensus.cl.df), row.names(consensus.cl.df))

color <- setNames(as.character(consensus.cl.df$cluster_color), row.names(consensus.cl.df))

pdf(file=paste0(output,"_dendrogram.pdf"))
dend.result <- build_dend(cl.med[,levels(cl.clean)], rank, color, nboot=30)

dend <- dend.result$dend

dend.labeled <- dend
labels(dend.labeled) <- consensus.cl.df[labels(dend), "cluster_label"]
plot(dend.labeled)
dev.off()
pdf(file=paste0(output,"_dendroheatmap.pdf"))
cl.cor <- dend.result$cl.cor
row.names(cl.cor) <- colnames(cl.cor) <- consensus.cl.df[row.names(cl.cor), "cluster_label"]
heatmap.3(cl.cor, Rowv=dend, Colv=dend, trace="none", col=heat.colors(100), cexRow=0.8, cexCol=0.8, breaks=c(-0.2, 0.2, seq(0.2, 1, length.out=99)))
dev.off()