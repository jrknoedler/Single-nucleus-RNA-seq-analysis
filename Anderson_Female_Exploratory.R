#!/usr/bin/env Rscript


output <- "Seurat/VMH_Anderson/Anderson_Exploratory_Female_Neuron"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)

options(future.globals.maxSize=1000*1024^2)

Female_Lactating_1.data <- Read10X(data.dir="andersonfemale/VMH_F/Lactating/2019_0820")
Female_Lactating_1 <- CreateSeuratObject(counts = Female_Lactating_1.data, project = "FemaleLactating1", min.cells=3, min.features=200)
Female_Lactating_1$sex <- "Female"
Female_Lactating_1$Hormone <- "Intact"
Female_Lactating_1$Condition <- "Lactating"

Female_Lactating_2.data <- Read10X(data.dir="andersonfemale/VMH_F/Lactating/2019_0924")
Female_Lactating_2 <- CreateSeuratObject(counts = Female_Lactating_2.data, project = "FemaleLactating2", min.cells=3, min.features=200)
Female_Lactating_2$sex <- "Female"
Female_Lactating_2$Hormone <- "Intact"
Female_Lactating_2$Condition <- "Lactating"

Female_Lactating_3.data <- Read10X(data.dir="andersonfemale/VMH_F/Lactating/2019_0925")
Female_Lactating_3 <- CreateSeuratObject(counts = Female_Lactating_3.data, project = "FemaleLactating3", min.cells=3, min.features=200)
Female_Lactating_3$sex <- "Female"
Female_Lactating_3$Hormone <- "Intact"
Female_Lactating_3$Condition <- "Lactating"

Female_Lactating_4.data <- Read10X(data.dir="andersonfemale/VMH_F/Lactating/2019_1126")
Female_Lactating_4 <- CreateSeuratObject(counts = Female_Lactating_4.data, project = "FemaleLactating4", min.cells=3, min.features=200)
Female_Lactating_4$sex <- "Female"
Female_Lactating_4$Hormone <- "Intact"
Female_Lactating_4$Condition <- "Lactating"

Female_Lactating_5.data <- Read10X(data.dir="andersonfemale/VMH_F/Lactating/2020_0117")
Female_Lactating_5 <- CreateSeuratObject(counts = Female_Lactating_5.data, project = "FemaleLactating5", min.cells=3, min.features=200)
Female_Lactating_5$sex <- "Female"
Female_Lactating_5$Hormone <- "Intact"
Female_Lactating_5$Condition <- "Lactating"

Female_Lactating_6.data <- Read10X(data.dir="andersonfemale/VMH_F/Lactating/2020_0326")
Female_Lactating_6 <- CreateSeuratObject(counts = Female_Lactating_6.data, project = "FemaleLactating6", min.cells=3, min.features=200)
Female_Lactating_6$sex <- "Female"
Female_Lactating_6$Hormone <- "Intact"
Female_Lactating_6$Condition <- "Lactating"


Female_Lactating_7.data <- Read10X(data.dir="andersonfemale/VMH_F/Lactating/2020_0405")
Female_Lactating_7 <- CreateSeuratObject(counts = Female_Lactating_7.data, project = "FemaleLactating7", min.cells=3, min.features=200)
Female_Lactating_7$sex <- "Female"
Female_Lactating_7$Hormone <- "Intact"
Female_Lactating_7$Condition <- "Lactating"


Female_Virgin_1.data <- Read10X(data.dir="andersonfemale/VMH_F/Virgin/2019_0822")
Female_Virgin_1 <- CreateSeuratObject(counts = Female_Virgin_1.data, project = "FemaleVirgin1", min.cells=3, min.features=200)
Female_Virgin_1$sex <- "Female"
Female_Virgin_1$Hormone <- "Intact"
Female_Virgin_1$Condition <- "Virgin"


Female_Virgin_2.data <- Read10X(data.dir="andersonfemale/VMH_F/Virgin/2019_0923")
Female_Virgin_2<- CreateSeuratObject(counts = Female_Virgin_2.data, project = "FemaleVirgin2", min.cells=3, min.features=200)
Female_Virgin_2$sex <- "Female"
Female_Virgin_2$Hormone <- "Intact"
Female_Virgin_2$Condition <- "Virgin"


Female_Virgin_3.data <- Read10X(data.dir="andersonfemale/VMH_F/Virgin/2019_1125")
Female_Virgin_3 <- CreateSeuratObject(counts = Female_Virgin_3.data, project = "FemaleVirgin3", min.cells=3, min.features=200)
Female_Virgin_3$sex <- "Female"
Female_Virgin_3$Hormone <- "Intact"
Female_Virgin_3$Condition <- "Virgin"

Female_Virgin_4.data <- Read10X(data.dir="andersonfemale/VMH_F/Virgin/2020_0118")
Female_Virgin_4 <- CreateSeuratObject(counts = Female_Virgin_4.data, project = "FemaleVirgin4", min.cells=3, min.features=200)
Female_Virgin_4$sex <- "Female"
Female_Virgin_4$Hormone <- "Intact"
Female_Virgin_4$Condition <- "Virgin"

Female_Virgin_5.data <- Read10X(data.dir="andersonfemale/VMH_F/Virgin/2020_0324")
Female_Virgin_5 <- CreateSeuratObject(counts = Female_Virgin_5.data, project = "FemaleVirgin5", min.cells=3, min.features=200)
Female_Virgin_5$sex <- "Female"
Female_Virgin_5$Hormone <- "Intact"
Female_Virgin_5$Condition <- "Virgin"

Female_Virgin_6.data <- Read10X(data.dir="andersonfemale/VMH_F/Virgin/2020_0325")
Female_Virgin_6 <- CreateSeuratObject(counts = Female_Virgin_6.data, project = "FemaleVirgin6", min.cells=3, min.features=200)
Female_Virgin_6$sex <- "Female"
Female_Virgin_6$Hormone <- "Intact"
Female_Virgin_6$Condition <- "Virgin"

Female_Virgin_7.data <- Read10X(data.dir="andersonfemale/VMH_F/Virgin/2020_0404")
Female_Virgin_7 <- CreateSeuratObject(counts = Female_Virgin_7.data, project = "FemaleVirgin7", min.cells=3, min.features=200)
Female_Virgin_7$sex <- "Female"
Female_Virgin_7$Hormone <- "Intact"
Female_Virgin_7$Condition <- "Virgin"

Female_Post_1.data <- Read10X(data.dir="andersonfemale/VMH_F/Post/2020_0209")
Female_Post_1 <- CreateSeuratObject(counts = Female_Post_1.data, project = "FemalePost1", min.cells=3, min.features=200)
Female_Post_1$sex <- "Female"
Female_Post_1$Hormone <- "Intact"
Female_Post_1$Condition <- "Post"

Female_Post_2.data <- Read10X(data.dir="andersonfemale/VMH_F/Post/2020_0304")
Female_Post_2 <- CreateSeuratObject(counts = Female_Post_2.data, project = "FemalePost2", min.cells=3, min.features=200)
Female_Post_2$sex <- "Female"
Female_Post_2$Hormone <- "Intact"
Female_Post_2$Condition <- "Post"



mySeurat <- merge(Female_Lactating_1, y=c(Female_Lactating_2,Female_Lactating_3,Female_Lactating_4,Female_Lactating_5,Female_Lactating_6,Female_Lactating_7,Female_Virgin_1,Female_Virgin_2,Female_Virgin_3,Female_Virgin_4,Female_Virgin_5,Female_Virgin_6,Female_Virgin_7,Female_Post_1,Female_Post_2), add.cell.ids=c("FemaleLactating1","FemaleLactating2","FemaleLactating3","FemaleLactating4","FemaleLactating5","FemaleLactating6","FemaleLactating7","FemaleVirgin1","FemaleVirgin2","FemaleVirgin3","FemaleVirgin4","FemaleVirgin5","FemaleVirgin6","FemaleVirgin7","FemalePost1","FemalePost2"), project="AndersonFemaleVMHMerged")
mySeurat
pdf(paste0(output,"RNAfeat.pdf"), height=7, width=28)
VlnPlot(mySeurat, features=("nCount_RNA"), split.by="orig.ident", pt.size=0)
dev.off()
pdf(paste0(output,"genenum.pdf"), height=7, width=28)
VlnPlot(mySeurat, features=c("nFeature_RNA"),split.by="orig.ident",pt.size=0)
dev.off()
#mySeurat[["percent.Malat1"]] <- PercentageFeatureSet(mySeurat, pattern = "Malat1")
mySeurat <- SCTransform(mySeurat, verbose=TRUE, vars.to.regress=c("orig.ident"))
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
mySeurat <- RunUMAP(mySeurat, reduction="pca", dim=1:50)
mySeurat <- FindNeighbors(mySeurat, reduction ="pca", dims=1:50)
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
mySeurat.markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.50, logfc.threshold = 0.50, test.use="MAST")
write.csv(mySeurat.markers, file=paste0(output,"_allposmarkers.csv"))
saveRDS(mySeurat, file=paste0(output, ".rds"))

