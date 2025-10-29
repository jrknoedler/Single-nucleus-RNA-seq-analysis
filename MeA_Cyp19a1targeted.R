#!/usr/bin/env Rscript

library(Seurat)

Primed <- readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered2_Primed_30pcs_Res1_Intact.rds")
Male <-  readRDS("Seurat/MeA_IndependentAnalysis/MeA_IndependentFiltered3_Intact_30pcs_1res_Intact.rds") 
Primed$sex <- "Primed"
Male$sex <- "Male" 
Primed <- subset(Primed, idents=c("16"))
Male <- subset(Male, idents=c("11"))
mySeurat <- merge(Primed, y=c(Male), add.cell.ids=c("Primed","Male"), project = "Merged") 
mySeurat$celltype.sex <- paste(mySeurat$sex) 
mySeurat$celltype <- "celltype.sex"
Idents(mySeurat) <- "celltype.sex"
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
markers <- FindMarkers(mySeurat, ident.1="Primed", ident.2="Male", test.use="MAST", min.pct=0, logfc.threshold=0)
write.csv(markers, file="Seurat/MeA_IndependentAnalysis/MeA_Neurogliaform_res1.5_MvP_SDEGs.csv")