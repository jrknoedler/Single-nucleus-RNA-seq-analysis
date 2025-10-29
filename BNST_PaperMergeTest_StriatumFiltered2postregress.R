#!/usr/bin/env Rscript

output <- "Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered_res1_sexclude_malat1regress_30pcs_Esr1filteredres1"
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_Striatumfiltered_res1_sexclude_malat1regress_30pcs_1.5res.rds")
mySeurat <- subset(mySeurat, idents=c("30","32"), invert=TRUE)
mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident","percent.Malat1"))
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
mySeurat <- RunPCA(mySeurat, features=hvg.final)
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:30)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:30)
mySeurat <- FindClusters(mySeurat, resolution=1)
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
pdf(paste0(output,"_UMAP_hormonelabel.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE, group.by="Hormone")
dev.off()
pdf(file=paste0(output,"Tac1_Vln.pdf"), width=40)
VlnPlot(mySeurat, features=c("Tac1"), pt.size=0, split.by="Hormone")
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
#top10 <- mySeurat.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
#DoHeatmap(mySeurat, features=top10$gene) + NoLegend()
#dev.off()
saveRDS(mySeurat, file=(paste0(output, ".rds")))
typemarkers <- read.table("DotPlotMarkerLists/BNST_Draft2_malat1regress.txt")
typemarkers <- unlist(typemarkers)
typemarkers <- as.matrix(typemarkers)
pdf(file=paste0(output, "Dotplot.pdf"), width=16)
DotPlot(mySeurat, features=typemarkers, dot.min=0.2, split.by="sex")
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