#!/usr/bin/env Rscript
rm(list=ls())


library(dplyr)
library(tidyverse)
library(plyr)

## Read list used for KB count
t2g <- read.table("/Volumes/SakuraSSD/ShahLab/Velocyte/kb_ref_vanilla2/t2g.txt", header = FALSE)
colnames(t2g) <- c("Transcript", "Gene ID", "Gene Name")
t2g <- t2g[-1]
t2g <- unique(t2g)


## VMH Pregnant

# Gene ID -> Gene Name 
spliced.genes <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/spliced/spliced.genes.txt", header = FALSE)
colnames(spliced.genes) <- "Gene ID"
data1 <- join(spliced.genes, t2g)
data1 <- data1[-1]
write.table(data1, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/Edit/spliced/spliced.genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

unspliced.genes <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/unspliced/unspliced.genes.txt", header = FALSE)
colnames(unspliced.genes) <- "Gene ID"
data2 <- join(unspliced.genes, t2g)
data2 <- data2[-1]
write.table(data2, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/Edit/unspliced/unspliced.genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# Rename cell barcode (set the same cell name using Naive Integration)
cell <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/spliced/spliced.barcodes.txt", header = FALSE)
cell$V1 <- sub("^","Pregnant_Pregnant_", cell$V1)
write.table(cell, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/Edit/spliced/spliced.barcodes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cell <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/unspliced/unspliced.barcodes.txt", header = FALSE)
cell$V1 <- sub("^","Pregnant_Pregnant_", cell$V1)
write.table(cell, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Preg_VMH/Edit/unspliced/unspliced.barcodes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


## VMH primed

# Gene ID -> Gene Name 
spliced.genes <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/spliced/spliced.genes.txt", header = FALSE)
colnames(spliced.genes) <- "Gene ID"
data1 <- join(spliced.genes, t2g)
data1 <- data1[-1]
write.table(data1, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/Edit/spliced/spliced.genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

unspliced.genes <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/unspliced/unspliced.genes.txt", header = FALSE)
colnames(unspliced.genes) <- "Gene ID"
data2 <- join(unspliced.genes, t2g)
data2 <- data2[-1]
write.table(data2, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/Edit/unspliced/unspliced.genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Rename cell barcode (set the same cell name using Naive Integration)
cell <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/spliced/spliced.barcodes.txt", header = FALSE)
cell$V1 <- sub("^","Ref_PrimedVMH_", cell$V1)
write.table(cell, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/Edit/spliced/spliced.barcodes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cell <- read.csv("/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/unspliced/unspliced.barcodes.txt", header = FALSE)
cell$V1 <- sub("^","Ref_PrimedVMH_", cell$V1)
write.table(cell, file = "/Volumes/SakuraSSD/ShahLab/Velocyte/Primed_VMH/Edit/unspliced/unspliced.barcodes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)