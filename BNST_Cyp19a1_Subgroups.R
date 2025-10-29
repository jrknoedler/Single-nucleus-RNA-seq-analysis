#!/usr/bin/env Rscript

library(Seurat)

Primed <- readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_Primed_30pcs_res1.2_lowthresh_Primed.rds")
Male <-  readRDS("Seurat/BNST_IndependentAnalysis/BNST_IndependentFiltered1_30pcs_res1_lowthresh_res1.2_Intact.rds") 
Primed$sex <- "Primed"
Male$sex <- "Male" 
Primed <- subset(Primed, idents=c("6"))
Male <- subset(Male, idents=c("12","13"))
mySeurat <- merge(Primed, y=c(Male), add.cell.ids=c("Primed","Male"), project = "Merged") 
mySeurat$celltype.sex <- paste(mySeurat$sex) 
mySeurat$celltype <- "celltype.sex"
Idents(mySeurat) <- "celltype.sex"
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
markers <- FindMarkers(mySeurat, ident.1="Primed", ident.2="Male", test.use="MAST", min.pct=0, logfc.threshold=0)
write.csv(markers, file="Seurat/BNST_IndependentAnalysis/BNST_Cyp19_Acvr1c_MvP_SDEGs.csv")