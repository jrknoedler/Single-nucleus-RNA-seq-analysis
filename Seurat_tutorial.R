#!/usr/bin/env Rscript


#load Seurat; packages below are generally optional but often useful
#For the Seurat functions, google them for list of all options
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viridis)
library(pheatmap)
library(patchwork)
# open the object
#
# Male cells = "Intact"; primed female cells ="Primed", 
mySeurat <- readRDS("path/to/file/filename.rds")

#Make sure you're using the RNA counts, not the weird normalized value used for the clustering
DefaultAssay(mySeurat) <- "RNA"

#long-normalize the RNA counts for better visualization if it isn't already log normalized
mySeurat <- NormalizeData(mySeurat)


#get some general information about your library
head(mySeurat[[]])

#The metadata includes the variables "sex" ("Male" or "Female") and "Hormone" ("Intact","Primed","Unprimed")

#make a violin plot for your favorite gene
VlnPlot(mySeurat, features=c("Tacr1"), pt.size=0)

#make a violin plot covered in ants
VlnPlot(mySeurat, features=c("Tacr1"), pt.size=1)

#make a violion plot of just your favorite gene in your favorite cluster, split by hormone status
VlnPlot(mySeurat, features=c("Tacr1"), idents=c("10"), split.by="Hormone", pt.size=0)

#make a stack of violin plots
VlnPlot(mySeurat, features=c("Tacr1","Esr1","Ar"), ncol=1, pt.size=0)

#make a UMAP
Dimplot(mySeurat, reduction="umap")

#Color your UMAP by sex
DimPlot(mySeurat, reduction="umap", group.by="sex")

#see where the cells expressing your gene are
FeaturePlot(mySeurat, features=c("Tacr1"), cols=c("gray","orange"))

#For any of the above, save the PDF for later (unless you're using Rstudio like a smart person and just export with the GUI
p <- VlnPlot(mySeurat, features=c("Tacr1","Esr1","Ar"), ncol=1, pt.size=0)
pdf(file="path/to/filename.pdf")
p
dev.off()



#figure out what gene are enriched in each cell type based on whatever criteria you like (may take 10-20 min)
markers <- FindAllMarkers(mySeurat, only.pos=TRUE, min.pct=0.25, logfc.threshold = 0.25)

#save it as a csv so you can actually look at it
write.csv(markers, file="path/to/output/filename.csv")

#Only look at male cells
maleSeurat <- subset(mySeurat, subset=sex=="Male")

#Filter a cluster you hate
filteredSeurat <- subset(mySeurat, idents=c("20"), invert=TRUE


#Rerun clustering if you think Joe did it wrong
DefaultAssay(mySeurat) <- "SCT"

mySeurat <- RunPCA(mySeurat)

#Use too many PCs (I used 30)

mySeurat <- FindNeighbors(mySeurat, dims=1:50)
#split the hell out of them (I used resolution 1.2)
mySeurat <- FindClusters(mySeurat, resolution=5)

mySeurat <- RunUMAP(mySeurat, dims=1:47)

#save your new oversplit Seurat object
saveRDS(mySeurat, file="path/to/filename.rds")