#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/VMH_IndependentAnalysis/VMH_Merge_Fig6storyboard_"

library(Seurat)
library(cowplot)
library(reticulate)
library(ggplot2)
library(dplyr)
library(MAST)
library(viper)
library(tidyverse) 
sessionInfo()
regdata <- file("ARACHNE/AllBulk/Expressionmtx.txt")
adjfile <- file("ARACHNE/AllBulk/network_3col.txt")
regul <- aracne2regulon(adjfile, regdata, format="3col", verbose=TRUE)

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_filtered3_Paperanalysis_30pcs.rds")
Male <- subset(x=mySeurat, subset=Hormone=="Intact")
Primed <- subset(x=mySeurat, subset=Hormone=="Primed")
Male.mat <- Male[["SCT"]]@scale.data
head(Male.mat)
#Male.mat <- as.matrix(Male.mat)
Primed.mat <- Primed[["SCT"]]@scale.data
#Primed.mat <- as.matrix(Primed.mat)
signatures <- viper::rowTtest(x=Male.mat, y=Primed.mat)
signatures <- (qnorm(signatures$p.value/2, lower.tail=FALSE) * sign(signatures$statistic))[,1]
null_model <- viper::ttestNull(x=Male.mat, y=Primed.mat, per=1, repos=F, verbose=T, cores=1, seed=1)