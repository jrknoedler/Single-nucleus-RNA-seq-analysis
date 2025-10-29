#!/usr/bin/env Rscript

output <- "Seurat/MeA_IndependentAnalysis/MeA_Fullmerge_redo_CCA"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
PrimedMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered2_Primed_30pcs_Res1.5_Intact.rds")
MaleMeA <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered3_Intact_30pcs_1.5res_Intact.rds")
UnprimedMeA <- readRDS("Seurat/UnprimedMeA_Exploratory/UnprimedMeA_emptyDrops_filtered1v3_30pcs_lowfdr_finalredo.rds")
#PrimedMeA <- SCTransform(PrimedMeA)
#MaleMeA <- SCTransform(MaleMeA)
#UnprimedMeA <- SCTransform(UnprimedMeA)
MaleMeA$sex <- "Male"
MaleMeA$Hormone <- "Intact"
PrimedMeA$sex <- "Female"
PrimedMeA$Hormone <- "Primed"
UnprimedMeA$sex <- "Female"
UnprimedMeA$Hormone <- "Unprimed"
list <- c(MaleMeA,UnprimedMeA,PrimedMeA)
MeA.features <- SelectIntegrationFeatures(object.list=list, nfeatures=3000)
MeA.list <- PrepSCTIntegration(list, anchor.features=MeA.features, verbose=TRUE)
MeA.anchors <- FindIntegrationAnchors(MeA.list, normalization.method="SCT", anchor.features=MeA.features, verbose=TRUE)
mySeurat <- IntegrateData(anchorset=MeA.anchors, normalization.method="SCT", verbose=FALSE)
genelist <- read.table("/home/groups/nirao/Bioinformatics3/Bioinformatics/ExcludeIDs.txt", header=FALSE)
genelist
genelist <- unlist(genelist)
genelist
genelist <- as.matrix(genelist)
genelist
mySeurat
head(mySeurat[[]])
hvg <- mySeurat@assays$SCT@var.features
hvg <- unlist(hvg)
hvg
hvg <- as.matrix(hvg)
hvg.final <- setdiff(hvg, genelist)
hvg.final
mySeurat <- RunPCA(mySeurat)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_UMAP_sexsplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_sexlabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="sex")
dev.off()
pdf(paste0(output,"_UMAP_hormonesplit.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, split.by="Hormone")
dev.off()
pdf(paste0(output,"_UMIvln.pdf"), width=12)
VlnPlot(mySeurat,features=c("nCount_RNA"), pt.size=1)
dev.off() 
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
pdf(paste0(output,"_Esr1vln.pdf"), width=16)
VlnPlot(mySeurat,features=c("Esr1"), pt.size=0)
dev.off() 
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=(paste0(output, ".rds")))
pdf(paste0(output,"_TopMarkerheatmap.pdf"))
top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))
pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.25)
dev.off()

Pct <- DotPlot(mySeurat, features=typemarkers)
sink(file=paste0(output,"sinkreally.txt"), append=FALSE, type=c("output","message"), split=FALSE)
Pct$data
sink()

data <- Pct$data
head(data)
pdf(file=paste0(output,"_rawdotplot.pdf"), width=16)
data %>% 
  filter(pct.exp > 25) %>% ggplot(aes(x=features.plot, y=factor(id, levels=rev(levels(factor(id)))), color=avg.exp, size=pct.exp)) + 
  geom_point() + 
  scale_color_viridis_c() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank())
dev.off()


# make data square to calculate euclidean distance
mat <- data %>% 
  select(-cluster) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = count) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$Gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix

pdf(file=paste0(output,"clusteredDotplot"))
dotplot <- gene_cluster %>% 
  mutate(Gene = factor(Gene, levels = clust$labels[clust$order])) %>% 
  #filter(pct.exp > 25) %>% 
  ggplot(aes(x=id, y=features.plot, color=avg.exp, size=pct.exp))  + 
  geom_point() + 
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')
dev.off()