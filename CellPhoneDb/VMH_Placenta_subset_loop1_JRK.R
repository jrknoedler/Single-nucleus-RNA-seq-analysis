#!/usr/bin/env Rscript
rm(list=ls())

library(Seurat)
library(dplyr)
library(tidyverse)
library(tibble)
library(tidyr)
library(rgeos)
library(reticulate)
library(future)
library(SeuratDisk)
library(scclusteval)
library(hdf5r)


###Load Brain + Placenta merged Seurat 

mySeurat <- readRDS("/scratch/users/tsakura/CellphoneDB/Integration_with_eLife_placenta/VMH2/VMH_Placenta_cell_type_HumanGenes_Add_ccell_type.rds")

 mySeurat <- SetIdent(mySeurat, value = mySeurat@meta.data$cell_type)



###specify some output destination
output <- "Subset_For_CellphoneDB/"

###Start your loop. Pay close attention to {} and ()

## Trophoblast

for (i in 0:31){  ###Create a vector that lists all clusters; if clusters are named (e.g. 'trophoblast', 'epithelial' you can also make a vector that lists these names
  try({    ###Embedding a 'try' loop makes the original loop keep going if one of the iterations fails

    sub1 <- subset(mySeurat, idents=c("Trophoblast")) 	
    sub2 <- subset(mySeurat, idents=c(i))
    sub <- merge(sub1, y=sub2)	
    
    save(sub, file=paste0(output, i, "_Trophoblast.rds")) ### variable name can be added to paste0 so each output file is unique
    
    SaveH5Seurat(sub, filename=paste0(output, i, "_Trophoblast_intermediate.h5Seurat"))
    Convert(paste0(output, i, "_Trophoblast_intermediate.h5Seurat"), dest = paste0(output, i, "_Trophoblast.h5ad"))
    
    sub@meta.data$cell_type  = sub@active.ident
    metadata <- sub@meta.data ##extract metadata
    metadata <- dplyr::select(metadata, cell_type)
    metadata$barcode_sample <- rownames(metadata)
    metadata = subset(metadata, select = c(barcode_sample, cell_type))
    metadata <- metadata[, c("barcode_sample", "cell_type")]
    write.table(metadata, file=paste0(output, i, "_Trophoblast_metadata.tsv"),
                row.names=FALSE, sep="\t")
  })
}



## Fetal Mesenchyme

for (i in 0:31){  ###Create a vector that lists all clusters; if clusters are named (e.g. 'trophoblast', 'epithelial' you can also make a vector that lists these names
  try({    ###Embedding a 'try' loop makes the original loop keep going if one of the iterations fails
    
    sub1 <- subset(mySeurat, idents=c("Fetal Mesenchyme")) 	
    sub2 <- subset(mySeurat, idents=c(i))
    sub <- merge(sub1, y=sub2)	
    
    save(sub, file=paste0(output, i, "_Fetal_Mesenchyme.rds")) ### variable name can be added to paste0 so each output file is unique
    
    SaveH5Seurat(sub, filename=paste0(output, i, "_Fetal_Mesenchyme_intermediate.h5Seurat"))
    Convert(paste0(output, i, "_Fetal_Mesenchyme_intermediate.h5Seurat"), dest = paste0(output, i, "_Fetal_Mesenchyme.h5ad"))
    
    sub@meta.data$cell_type  = sub@active.ident
    metadata <- sub@meta.data ##extract metadata
    metadata <- dplyr::select(metadata, cell_type)
    metadata$barcode_sample <- rownames(metadata)
    metadata = subset(metadata, select = c(barcode_sample, cell_type))
    metadata <- metadata[, c("barcode_sample", "cell_type")]
    write.table(metadata, file=paste0(output, i, "_Fetal_Mesenchyme_metadata.tsv"),
                row.names=FALSE, sep="\t")
  })
}


## Endothelial

for (i in 0:31){  ###Create a vector that lists all clusters; if clusters are named (e.g. 'trophoblast', 'epithelial' you can also make a vector that lists these names
  try({    ###Embedding a 'try' loop makes the original loop keep going if one of the iterations fails
    
    sub1 <- subset(mySeurat, idents=c("Endothelial")) 	
    sub2 <- subset(mySeurat, idents=c(i))
    sub <- merge(sub1, y=sub2)
    
    save(sub, file=paste0(output, i, "_Endothelial.rds")) ### variable name can be added to paste0 so each output file is unique
    
    SaveH5Seurat(sub, filename=paste0(output, i, "_Endothelial_intermediate.h5Seurat"))
    Convert(paste0(output, i, "_Endothelial_intermediate.h5Seurat"), dest = paste0(output, i, "_Endothelial.h5ad"))
    
    sub@meta.data$cell_type  = sub@active.ident
    metadata <- sub@meta.data ##extract metadata
    metadata <- dplyr::select(metadata, cell_type)
    metadata$barcode_sample <- rownames(metadata)
    metadata = subset(metadata, select = c(barcode_sample, cell_type))
    metadata <- metadata[, c("barcode_sample", "cell_type")]
    write.table(metadata, file=paste0(output, i, "_Endothelial_metadata.tsv"),
                row.names=FALSE, sep="\t")
  })
}

## Blood Cells

for (i in 0:31){  ###Create a vector that lists all clusters; if clusters are named (e.g. 'trophoblast', 'epithelial' you can also make a vector that lists these names
  try({    ###Embedding a 'try' loop makes the original loop keep going if one of the iterations fails
    
    sub1 <- subset(mySeurat, idents=c("Blood Cells")) 	
    sub2 <- subset(mySeurat, idents=c(i))
    sub <- merge(sub1, y=sub2)

    save(sub, file=paste0(output, i, "_Blood_Cells.rds")) ### variable name can be added to paste0 so each output file is unique
    
    SaveH5Seurat(sub, filename=paste0(output, i, "_Blood_Cells_intermediate.h5Seurat"))
    Convert(paste0(output, i, "_Blood_Cells_intermediate.h5Seurat"), dest = paste0(output, i, "_Blood_Cells.h5ad"))
    
    sub@meta.data$cell_type  = sub@active.ident
    metadata <- sub@meta.data ##extract metadata
    metadata <- dplyr::select(metadata, cell_type)
    metadata$barcode_sample <- rownames(metadata)
    metadata = subset(metadata, select = c(barcode_sample, cell_type))
    metadata <- metadata[, c("barcode_sample", "cell_type")]
    write.table(metadata, file=paste0(output, i, "_Blood_Cells_metadata.tsv"),
                row.names=FALSE, sep="\t")
  })
}

## Decidual Stroma

for (i in 0:31){  ###Create a vector that lists all clusters; if clusters are named (e.g. 'trophoblast', 'epithelial' you can also make a vector that lists these names
  try({    ###Embedding a 'try' loop makes the original loop keep going if one of the iterations fails
    
    sub1 <- subset(mySeurat, idents=c("Decidual Stroma")) 	
    sub2 <- subset(mySeurat, idents=c(i))
    sub <- merge(sub1, y=sub2)

    save(sub, file=paste0(output, i, "_Decidual_Stroma.rds")) ### variable name can be added to paste0 so each output file is unique
    
    SaveH5Seurat(sub, filename=paste0(output, i, "_Decidual_Stroma_intermediate.h5Seurat"))
    Convert(paste0(output, i, "_Decidual_Stroma_intermediate.h5Seurat"), dest = paste0(output, i, "_Decidual_Stroma.h5ad"))
    
    sub@meta.data$cell_type  = sub@active.ident
    metadata <- sub@meta.data ##extract metadata
    metadata <- dplyr::select(metadata, cell_type)
    metadata$barcode_sample <- rownames(metadata)
    metadata = subset(metadata, select = c(barcode_sample, cell_type))
    metadata <- metadata[, c("barcode_sample", "cell_type")]
    write.table(metadata, file=paste0(output, i, "_Decidual_Stroma_metadata.tsv"),
                row.names=FALSE, sep="\t")
  })
}
