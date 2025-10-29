#!/usr/bin/env Rscript
library(Seurat)
library(rgeos)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(future)
library(SeuratDisk)
library(scclusteval)
library(hdf5r)


data <- readRDS("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH2/VMH_Placenta_cell_type_merged_NormalizedRNA.rds")

SaveH5Seurat(data, filename="/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH2/VMH_Placenta_cell_type_merged_NormalizedRNA_intermediate.h5Seurat")

Convert("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH2/VMH_Placenta_cell_type_merged_NormalizedRNA_intermediate.h5Seurat", dest = "/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH2/VMH_Placenta_cell_type_merged_NormalizedRNA.h5ad")