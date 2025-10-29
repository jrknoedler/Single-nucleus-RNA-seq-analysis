#!/usr/bin/env Rscript



output <- "LIGERoutput/VMH/VMH_Liger_mk1"
library(liger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(cowplot)

Primed <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1_Primed_30pcs_res1_Primed.rds")
Primed <- RenameCells(Primed, add.cell.id="Primed")
Male <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered1_Intact_30pcs_res1_Primed.rds")
Male <- RenameCells(Male, add.cell.id="Male")
Unprimed <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_IndependentFiltered3_Unprimed_30pcs_res1.rds")
Unprimed <- RenameCells(Unprimed, add.cell.id="Unprimed")
Male$sex <- "Male"
Male$Hormone <- "Intact"
Primed$sex <- "Female"
Primed$Hormone <- "Primed"
Unprimed$sex <- "Female"
Unprimed$Hormone <- "Unprimed"

Primed_Mat <- Primed@assays[["RNA"]]@counts
Primed_Mat <- as.matrix(Primed_Mat)
Male_Mat <- Male@assays[["RNA"]]@counts
Male_Mat <- as.matrix(Male_Mat)
Unprimed_Mat <- Unprimed@assays[["RNA"]]@counts
Unprimed_Mat <- as.matrix(Unprimed_Mat)

liger <- createLiger(list(Male=Male_Mat, Primed=Primed_Mat, Unprimed=Unprimed_Mat))

liger <- normalize(liger)
liger <- selectGenes(liger)
liger <- scaleNotCenter(liger)
liger <- optimizeALS(liger, k=40, lambda=1)

liger <- quantile_norm(liger)

liger <- runUMAP(liger, distance='cosine', n_neighbors=30, min_dist=0.3)

pdf(file=paste0(output,"_umap.pdf"))
all.plots <- plotByDatasetAndCluster(liger, axis.labels=c('UMAP 1', 'UMAP 2'), return.plots=T)
all.plots[[1]] + all.plots[[2]]
dev.off()
pdf(file=paste0(output,"_factor1loadings.pdf"))
gene_loadings <- plotGeneLoadings(liger, num.genes = 10, do.spec.plot = F, return.plots=TRUE)
gene_loadings[[1]]
dev.off()
pdf(file=paste0(output,"_factor2loadings.pdf"))
gene_loadings <- plotGeneLoadings(liger, num.genes = 10, do.spec.plot = F,return.plots=TRUE)
gene_loadings[[2]]
dev.off()
pdf(file=paste0(output,"_factor3loadings.pdf"))
gene_loadings <- plotGeneLoadings(liger, num.genes = 10, do.spec.plot = F,return.plots=TRUE)
gene_loadings[[3]]
dev.off()
cluster.results <- runWilcoxon(liger, compare.method="clusters")
head(cluster.results)
datasets.results <- runWilcoxon(liger, compare.method="datasets")
head(datasets.results)
#pdf(file=paste0(output,"_riverplot.pdf"))
#makeRiverplot(liger, Male, Primed, min.frac=0.025)
#dev.off()
cluster.results <- cluster.results[cluster.results$padj < 0.05,]
datasets.results <- datasets.results[datasets.results$padj < 0.05,]
write.csv(cluster.results, file=paste0(output,"_clusterresults.csv"))
write.csv(datasets.results, file=paste0(output,"_datasetsresults.csv"))
saveRDS(liger, file=paste0(output,"_ligertrial.rds"))