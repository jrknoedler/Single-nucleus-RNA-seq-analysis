library(Seurat)

output <- 

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/

data <- mySeurat@assays[["RNA"]]@counts
clusters <- mySeurat@Idents

Downsample <- SampleUMI(data, max.umi=3000, upsample=FALSE, verbose=TRUE)

Seurat.ds <- CreateSeuratObject(Downsample)
