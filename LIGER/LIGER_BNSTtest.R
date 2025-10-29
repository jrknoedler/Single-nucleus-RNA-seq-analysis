#!/usr/bin/env Rscript



output <- "LIGERoutput/BNST/BNST_Liger_100k"
library(hdf5r)
library(rliger)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(cowplot)

Primed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_Primed_40pcs_res1.5_lowthresh_Primed.rds")
Primed <- RenameCells(Primed, add.cell.id="Primed")
Male <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_40pcs_res1_lowthresh_res1.5_Intact.rds")
Male <- RenameCells(Male, add.cell.id="Male")
Unprimed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered2_Unprimed_30pcs_nothresh_res1__Unprimed.rds")
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
liger <- optimizeALS(liger, k=100, lambda=1)

liger <- quantile_norm(liger, eps=0.5)

liger <- runUMAP(liger, distance='cosine', n_neighbors=30, min_dist=0.3)

pdf(file=paste0(output,"_umap.pdf"))
all.plots <- plotByDatasetAndCluster(liger, axis.labels=c('UMAP 1', 'UMAP 2'), return.plots=T)
all.plots[[1]] + all.plots[[2]]
dev.off()
pdf(file=paste0(output,"_geneloadings.pdf"))
plotGeneLoadings(liger, return.plots=TRUE)
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
