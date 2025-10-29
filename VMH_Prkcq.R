#!/usr/bin/env Rscript

library(Seurat)

Primed <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Primed_CuratedPrimed_Merged.rds")
Male <- readRDS("Seurat/VMH_IndependentAnalysis/VMH_Intact_CuratedIntact_Merged.rds") 
Primed$sex <- "Primed"
Male$sex <- "Male" 
Primed <- subset(Primed, idents=c("Prkcq"))
Male <- subset(Male, idents=c("Prkcq"))
mySeurat <- merge(Primed, y=c(Male), add.cell.ids=c("Primed","Male"), project = "Merged") 
mySeurat$celltype.sex <- paste(mySeurat$sex) 
mySeurat$celltype <- "celltype.sex"
Idents(mySeurat) <- "celltype.sex"
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
markers <- FindMarkers(mySeurat, ident.1="Primed", ident.2="Male", test.use="MAST", min.pct=0, logfc.threshold=0)
write.csv(markers, file="Seurat/VMH_IndependentAnalysis/VMH_Prkcq_SDEGs.csv")