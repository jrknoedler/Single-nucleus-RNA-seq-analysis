library(Seurat)

mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_30pcs_res1.2_malatregressexcludeFINAL.rds")

idents <- mySeurat@active.ident
hormone <- mySeurat$Hormone

umap <- mySeurat[["umap"]]@cell.embeddings
umap <- as.data.frame(umap)

final <- cbind(umap, idents, hormone)
head(final)

write.csv(final, file="ForChungha.csv")
