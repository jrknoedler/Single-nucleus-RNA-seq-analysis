#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "GEO_Upload_2/POA_Unfiltered_"

library(BUSpaRse)
library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
Primed.data <- read_count_output("Kallisto_Bus_Output/PrimedPOA/combined_counts_em/", name="output", tcc=FALSE)
Primed <- CreateSeuratObject(counts = Primed.data, project = "PrimedFemale", min.cells=3, min.features=500)
Primedcounts <- Primed@assays$RNA@counts
write.table(Primedcounts, file=paste0(output,"Primed.txt"))
Intact.data <- read_count_output("Kallisto_Bus_Output/IntactPOA/combined_counts_em/", name="output", tcc=FALSE)
Intact <- CreateSeuratObject(counts = Intact.data, project = "IntactMale", min.cells=3, min.features=500)
Intactcounts <- Intact@assays$RNA@counts
write.table(Intactcounts, file=paste0(output,"Intact.txt"))
Unprimed.data <- read_count_output("Kallisto_Bus_Output/UnprimedPOA/combined_counts_em/", name="output", tcc=FALSE)
Unprimed <- CreateSeuratObject(counts = Unprimed.data, project = "UnprimedFemale", min.cells=3, min.features=500)
Unprimedcounts <- Unprimed@assays$RNA@counts
write.table(Unprimedcounts, file=paste0(output,"Unprimed.txt"))

