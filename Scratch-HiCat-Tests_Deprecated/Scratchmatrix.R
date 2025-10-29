#!/usr/bin/env Rscript



output <- "Scrattch_hicat/BNST/BNSTMerged_scrattch_hicat_May"

mySeurat <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_Fullmerge_Test_prepaperreclusterRound1_sexclude_malat1exclude_malat1regress.rds")
cells <- rownames(mySeurat[[]])
data <- mySeurat@assays[["RNA"]]@counts
data <- as.matrix(data)
saveRDS(data, file=paste0(output,"BNSTMatrix.rds"))