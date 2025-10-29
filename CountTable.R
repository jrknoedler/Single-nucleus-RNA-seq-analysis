library(Seurat)

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_filtered3_Paperanalysis_20pcs.rds")

mySeurat <- NormalizeData(mySeurat)
data <- mySeurat@assays[["RNA"]]@counts
write.table(data, file="ARACHNE/VMH/VMH_Matrixtest.txt")