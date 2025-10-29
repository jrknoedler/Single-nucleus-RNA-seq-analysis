#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_Anderson/Anderson_Exploratory"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)

Female_Control_1.data <- Read10X(data.dir="Anderson_Data/Female_Control_1")
Female_Control_1 <- CreateSeuratObject(counts = Female_Control_1.data, project = "FemaleControl1", min.cells=3, min.features=200)
Female_Control_1$sex <- "Female"
Female_Control_1$Hormone <- "Intact"
Female_Control_1$Condition <- "Control"

Female_Fnr_1.data <- Read10X(data.dir="Anderson_Data/Female_Fnr_1")
Female_Fnr_1 <- CreateSeuratObject(counts = Female_Fnr_1.data, project = "FemaleFnr1", min.cells=3, min.features=200)
Female_Fnr_1$sex <- "Female"
Female_Fnr_1$Hormone <- "Intact"
Female_Fnr_1$Condition <- "Non-receptive"

Female_Fnr_2.data <- Read10X(data.dir="Anderson_Data/Female_Fnr_2")
Female_Fnr_2 <- CreateSeuratObject(counts = Female_Fnr_2.data, project = "FemaleFnr2", min.cells=3, min.features=200)
Female_Fnr_2$sex <- "Female"
Female_Fnr_2$Hormone <- "Intact"
Female_Fnr_2$Condition <- "Non-receptive"

Female_Fnr_3.data <- Read10X(data.dir="Anderson_Data/Female_Fnr_3")
Female_Fnr_3 <- CreateSeuratObject(counts = Female_Fnr_3.data, project = "FemaleFnr3", min.cells=3, min.features=200)
Female_Fnr_3$sex <- "Female"
Female_Fnr_3$Hormone <- "Intact"
Female_Fnr_3$Condition <- "Non-receptive"

Female_Plain_1.data <- Read10X(data.dir="Anderson_Data/Female_Plain_1")
Female_Plain_1 <- CreateSeuratObject(counts = Female_Plain_1.data, project = "FemalePlain1", min.cells=3, min.features=200)
Female_Fnr_3$sex <- "Female"
Female_Fnr_3$Hormone <- "Intact"
Female_Fnr_3$Condition <- "Plain"

Male_Aggression_1.data <- Read10X(data.dir="Anderson_Data/Male_Aggression_1")
Male_Aggression_1 <- CreateSeuratObject(counts = Male_Aggression_1.data, project = "MaleAggression1", min.cells=3, min.features=200)
Male_Aggression_1$sex <- "Male"
Male_Aggression_1$Hormone <- "Intact"
Male_Aggression_1$Condition <- "Aggression"


Male_Aggression_2.data <- Read10X(data.dir="Anderson_Data/Male_Aggression_2")
Male_Aggression_2 <- CreateSeuratObject(counts = Male_Aggression_2.data, project = "MaleAggression2", min.cells=3, min.features=200)
Male_Aggression_2$sex <- "Male"
Male_Aggression_2$Hormone <- "Intact"
Male_Aggression_2$Condition <- "Aggression"


Male_Aggression_3.data <- Read10X(data.dir="Anderson_Data/Male_Aggression_3")
Male_Aggression_3 <- CreateSeuratObject(counts = Male_Aggression_3.data, project = "MaleAggression3", min.cells=3, min.features=200)
Male_Aggression_3$sex <- "Male"
Male_Aggression_3$Hormone <- "Intact"
Male_Aggression_3$Condition <- "Aggression"


Male_Control_1.data <- Read10X(data.dir="Anderson_Data/Male_Control_1")
Male_Control_1 <- CreateSeuratObject(counts = Male_Control_1.data, project = "MaleControl1", min.cells=3, min.features=200)
Male_Control_1$sex <- "Male"
Male_Control_1$Hormone <- "Intact"
Male_Control_1$Condition <- "Control"


Male_Control_2.data <- Read10X(data.dir="Anderson_Data/Male_Control_2")
Male_Control_2 <- CreateSeuratObject(counts = Male_Control_2.data, project = "MaleControl2", min.cells=3, min.features=200)
Male_Control_2$sex <- "Male"
Male_Control_2$Hormone <- "Intact"
Male_Control_2$Condition <- "Control"

Male_Control_3.data <- Read10X(data.dir="Anderson_Data/Male_Control_3")
Male_Control_3 <- CreateSeuratObject(counts = Male_Control_3.data, project = "FemaleControl1", min.cells=3, min.features=200)
Male_Control_3$sex <- "Male"
Male_Control_3$Hormone <- "Intact"
Male_Control_3$Condition <- "Control"

Male_Mating_1.data <- Read10X(data.dir="Anderson_Data/Male_Mating_1")
Male_Mating_1 <- CreateSeuratObject(counts = Male_Mating_1.data, project = "MaleMating1", min.cells=3, min.features=200)
Male_Mating_1$sex <- "Male"
Male_Mating_1$Hormone <- "Intact"
Male_Mating_1$Condition <- "Control"

Male_MF_CI_dangled_1.data <- Read10X(data.dir="Anderson_Data/Male_MF_CI_dangled_1")
Male_MF_CI_dangled_1 <- CreateSeuratObject(counts = Male_MF_CI_dangled_1.data, project = "MaleMFCIdangled1", min.cells=3, min.features=200)
Male_MF_CI_dangled_1$sex <- "Male"
Male_MF_CI_dangled_1$Hormone <- "Intact"
Male_MF_CI_dangled_1$Condition <- "Male MF CI Dangled"

Male_MF_CI_pencilcup.data <- Read10X(data.dir="Anderson_Data/Male_MF_CI_pencilcup")
Male_MF_CI_pencilcup <- CreateSeuratObject(counts = Male_MF_CI_pencilcup.data, project = "MaleMFCIpencilcup1", min.cells=3, min.features=200)
Male_MF_CI_pencilcup$sex <- "Male"
Male_MF_CI_pencilcup$Hormone <- "Intact"
Male_MF_CI_pencilcup$Condition <- "Male MF CI Pencil"

Male_MM_CI_dangled_1.data <- Read10X(data.dir="Anderson_Data/Male_MM_CI_dangled_1")
Male_MM_CI_dangled_1 <- CreateSeuratObject(counts = Male_MM_CI_dangled_1.data, project = "MaleMMCIdangled1", min.cells=3, min.features=200)
Male_MM_CI_dangled_1$sex <- "Male"
Male_MM_CI_dangled_1$Hormone <- "Intact"
Male_MM_CI_dangled_1$Condition <- "Male MM CI Dangled"

Male_MM_CI_dangled_2.data <- Read10X(data.dir="Anderson_Data/Male_MM_CI_dangled_2")
Male_MM_CI_dangled_2 <- CreateSeuratObject(counts = Male_MM_CI_dangled_2.data, project = "MaleMMCIdangled2", min.cells=3, min.features=200)
Male_MM_CI_dangled_2$sex <- "Male"
Male_MM_CI_dangled_2$Hormone <- "Intact"
Male_MM_CI_dangled_2$Condition <- "Male MM CI Dangled"

Male_MM_CI_pencilcup.data <- Read10X(data.dir="Anderson_Data/Male_MM_CI_pencilcup")
Male_MM_CI_pencilcup <- CreateSeuratObject(counts = Male_MM_CI_pencilcup.data, project = "MaleMMCIpencilcup", min.cells=3, min.features=200)
Male_MM_CI_pencilcup$sex <- "Male"
Male_MM_CI_pencilcup$Hormone <- "Intact"
Male_MM_CI_pencilcup$Condition <- "Male MM CI Pencilcup"

Male_Other_1.data <- Read10X(data.dir="Anderson_Data/Male_Other_1")
Male_Other_1 <- CreateSeuratObject(counts = Male_Other_1.data, project = "MaleOther1", min.cells=3, min.features=200)
Male_Other_1$sex <- "Male"
Male_Other_1$Hormone <- "Intact"
Male_Other_1$Condition <- "Male Other"

Male_Other_2.data <- Read10X(data.dir="Anderson_Data/Male_Other_2")
Male_Other_2 <- CreateSeuratObject(counts = Male_Other_2.data, project = "MaleOther2", min.cells=3, min.features=200)
Male_Other_2$sex <- "Male"
Male_Other_2$Hormone <- "Intact"
Male_Other_2$Condition <- "Male Other"

Male_Other_3.data <- Read10X(data.dir="Anderson_Data/Male_Other_3")
Male_Other_3 <- CreateSeuratObject(counts = Male_Other_3.data, project = "MaleOther3", min.cells=3, min.features=200)
Male_Other_3$sex <- "Male"
Male_Other_3$Hormone <- "Intact"
Male_Other_3$Condition <- "Male Other"

Male_Other_4.data <- Read10X(data.dir="Anderson_Data/Male_Other_4")
Male_Other_4 <- CreateSeuratObject(counts = Male_Other_4.data, project = "MaleOther4", min.cells=3, min.features=200)
Male_Other_4$sex <- "Male"
Male_Other_4$Hormone <- "Intact"
Male_Other_4$Condition <- "Male Other"

Male_Plain_1.data <- Read10X(data.dir="Anderson_Data/Male_Plain_1")
Male_Plain_1 <- CreateSeuratObject(counts = Male_Plain_1.data, project = "MalePlain1", min.cells=3, min.features=200)
Male_Plain_1$sex <- "Male"
Male_Plain_1$Hormone <- "Intact"
Male_Plain_1$Condition <- "Male Plain"

Male_Plain_2.data <- Read10X(data.dir="Anderson_Data/Male_Plain_2")
Male_Plain_2 <- CreateSeuratObject(counts =Male_Plain_2.data, project = "MalePlain2", min.cells=3, min.features=200)
Male_Plain_2$sex <- "Male"
Male_Plain_2$Hormone <- "Intact"
Male_Plain_2$Condition <- "Male Plain"

Male_Social_Fear_Group.data <- Read10X(data.dir="Anderson_Data/Male_Social_Fear_Group")
Male_Social_Fear_Group <- CreateSeuratObject(counts = Male_Social_Fear_Group.data, project = "MaleSocFearGroup", min.cells=3, min.features=200)
Male_Social_Fear_Group$sex <- "Male"
Male_Social_Fear_Group$Hormone <- "Intact"
Male_Social_Fear_Group$Condition <- "Male Social Fear Group"

Male_Social_Fear_Single_1.data <- Read10X(data.dir="Anderson_Data/Male_Social_Fear_Single_1")
Male_Social_Fear_Single_1 <- CreateSeuratObject(counts = Male_Social_Fear_Single_1.data, project = "MaleSocFearSingle1", min.cells=3, min.features=200)
Male_Social_Fear_Single_1$sex <- "Male"
Male_Social_Fear_Single_1$Hormone <- "Intact"
Male_Social_Fear_Single_1$Condition <- "Male Social Fear Single"

Male_Social_Fear_Single_2.data <- Read10X(data.dir="Anderson_Data/Male_Social_Fear_Single_2")
Male_Social_Fear_Single_2 <- CreateSeuratObject(counts = Male_Social_Fear_Single_2.data, project = "MaleSocFearSingle2", min.cells=3, min.features=200)
Male_Social_Fear_Single_2$sex <- "Male"
Male_Social_Fear_Single_2$Hormone <- "Intact"
Male_Social_Fear_Single_2$Condition <- "Male Social Fear Single"

mySeurat <- merge(Female_Control_1, y=c(Female_Fnr_1,Female_Fnr_2,Female_Fnr_3,Female_Plain_1,Male_Aggression_1,Male_Aggression_2,Male_Aggression_3,Male_Control_1,Male_Control_2,Male_Control_3,Male_Mating_1,Male_MF_CI_dangled_1,Male_MF_CI_pencilcup,Male_MM_CI_dangled_1,Male_MM_CI_dangled_2,Male_MM_CI_pencilcup,Male_Other_1,Male_Other_2,Male_Other_3,Male_Other_4,Male_Plain_1,Male_Plain_2,Male_Social_Fear_Group,Male_Social_Fear_Single_1,Male_Social_Fear_Single_2), add.cell.ids=c("Female_Control_1","Female_Fnr_1","Female_Fnr_2","Female_Fnr_3","Female_Plain_1","Male_Aggression_1","Male_Aggression_2","Male_Aggression_3","Male_Control_1","Male_Control_2","Male_Control_3","Male_Mating_1","Male_MF_CI_dangled_1","Male_MF_CI_pencilcup","Male_MM_CI_dangled_1","Male_MM_CI_dangled_2","Male_MM_CI_pencilcup","Male_Other_1","Male_Other_2","Male_Other_3","Male_Other_4","Male_Plain_1","Male_Plain_2","Male_Social_Fear_Group","Male_Social_Fear_Single_1","Male_Social_Fear_Single_2"), project="AndersonVMHMerged")
mySeurat
pdf(paste0(output,"RNAfeat.pdf"))
VlnPlot(mySeurat, features=("nCount_RNA"), split.by="Condition", pt.size=0)
dev.off()
pdf(paste0(output,"genenum.pdf"))
VlnPlot(mySeurat, features=c("nFeature_RNA"),split.by="Condition",pt.size=0)
dev.off()
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
mySeurat <- FindClusters(mySeurat, resolution=1.2)
pdf(paste0(output,"_UMAP.pdf"))
DimPlot(mySeurat, reduction="umap", label=TRUE)
dev.off()
pdf(paste0(output,"_MarkerVln.pdf"), width=40, height=60)
VlnPlot(mySeurat, features=c("Gad1","Slc17a6","Esr1","Cckar","Mbp","Plp1","Cxcr1","Inppd5"), ncol=1, pt.size=0)
dev.off()
SeuratScaled <- ScaleData(mySeurat)
pdf(paste0(output,"_Markerdot.pdf"), width=40, height=40)
DotPlot(SeuratScaled, features=c("Gad1","Slc17a6","Esr1","Cckar","Mbp","Plp1","Cxcr1","Inppd5"))
dev.off()
DefaultAssay(mySeurat) <- "RNA"
mySeurat <- NormalizeData(mySeurat)
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=paste0(output, ".rds"))

