
library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)



Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedPOAattempt/combined_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=1, min.features=200)
Primed

pdf(file="PrimedPOA_Kiss1Search_Gm28040.pdf")
VlnPlot(Primed, features=c("Gm28040"))
dev.off()
pdf(file="PrimedPOA_Kiss1Search_Golt1a.pdf")
VlnPlot(Primed, features=c("Golt1a"))
dev.off()
pdf(file="PrimedPOA_Kiss1Search__Kiss1.pdf")
VlnPlot(Primed, features=c("Kiss1"))
dev.off()