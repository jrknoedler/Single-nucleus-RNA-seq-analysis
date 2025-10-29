#!/usr/bin/env Rscript

library(Seurat)

Primed <- readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentFiltered2_30pcs_res1_Primed.rds")
Male <-  readRDS("Seurat/POA_IndependentAnalysis/POA_IndependentFiltered7_30pcsredo_res1_Intact.rds") 
Primed$sex <- "Primed"
Male$sex <- "Male" 
Primed <- subset(Primed, idents=c("3","2","4","16","10","21","24","20"))
Male <- subset(Male, idents=c("4","16","18","15","6","17","21","9","13"))
mySeurat <- merge(Primed, y=c(Male), add.cell.ids=c("Primed","Male"), project = "Merged") 
mySeurat$celltype.sex <- paste(mySeurat$sex) 
mySeurat$celltype <- "celltype.sex"
Idents(mySeurat) <- "celltype.sex"
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
markers <- FindMarkers(mySeurat, ident.1="Primed", ident.2="Male", test.use="MAST", min.pct=0, logfc.threshold=0)
write.csv(markers, file="Seurat/POA_IndependentAnalysis/POA_Supergroup1_MvP_SDEGs.csv")