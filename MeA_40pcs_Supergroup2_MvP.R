#!/usr/bin/env Rscript

library(Seurat)

Primed <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered2_Primed_40pcs_Res1_Intact.rds")
Male <-  readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered3_Intact_40pcs_1res_Intact.rds") 
Primed$sex <- "Primed"
Male$sex <- "Male" 
Primed <- subset(Primed, idents=c("16","4","5"))
Male <- subset(Male, idents=c("5","3","22"))
mySeurat <- merge(Primed, y=c(Male), add.cell.ids=c("Primed","Male"), project = "Merged") 
mySeurat$celltype.sex <- paste(mySeurat$sex) 
mySeurat$celltype <- "celltype.sex"
Idents(mySeurat) <- "celltype.sex"
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
markers <- FindMarkers(mySeurat, ident.1="Primed", ident.2="Male", min.pct=0, logfc.threshold=0, test.use="MAST")
write.csv(markers, file="Seurat/MeA_IndependentAnalysis/MeA_Supergroup2_minimal_res1_40pcs_MvP_SDEGsdefault.csv")