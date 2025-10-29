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

##### The three lines of code between the comments are sufficient to convert
data <- readRDS("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset2_VMH_Placenta_HumanGenes.rds")

SaveH5Seurat(data, filename="/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset2_VMH_Placenta_HumanGenes_intermediate.h5Seurat")

Convert("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset2_VMH_Placenta_HumanGenes_intermediate.h5Seurat", dest = "/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset2_VMH_Placenta_HumanGenes.h5ad")
#####
data <- readRDS("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset3_VMH_Placenta_HumanGenes.rds")

SaveH5Seurat(data, filename="/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset3_VMH_Placenta_HumanGenes_intermediate.h5Seurat")

Convert("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset3_VMH_Placenta_HumanGenes_intermediate.h5Seurat", dest = "/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset3_VMH_Placenta_HumanGenes.h5ad")

data <- readRDS("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset4_VMH_Placenta_HumanGenes.rds")

SaveH5Seurat(data, filename="/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset4_VMH_Placenta_HumanGenes_intermediate.h5Seurat")

Convert("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset4_VMH_Placenta_HumanGenes_intermediate.h5Seurat", dest = "/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset4_VMH_Placenta_HumanGenes.h5ad")

data <- readRDS("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset5_VMH_Placenta_HumanGenes.rds")

SaveH5Seurat(data, filename="/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset5_VMH_Placenta_HumanGenes_intermediate.h5Seurat")

Convert("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset5_VMH_Placenta_HumanGenes_intermediate.h5Seurat", dest = "/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset5_VMH_Placenta_HumanGenes.h5ad")

data <- readRDS("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset6_VMH_Placenta_HumanGenes.rds")

SaveH5Seurat(data, filename="/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset6_VMH_Placenta_HumanGenes_intermediate.h5Seurat")

Convert("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset6_VMH_Placenta_HumanGenes_intermediate.h5Seurat", dest = "/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH_separate1/subset6_VMH_Placenta_HumanGenes.h5ad")