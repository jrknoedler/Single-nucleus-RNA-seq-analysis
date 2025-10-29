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

regdata <- file("ARACHNE/AllBulkvst/Expressionmtx.txt")
adjfile <- file("ARACHNE/AllBulkvst/network_3col.txt")
regul <- aracne2regulon(adjfile, regdata, format="3col", verbose=TRUE)

mySeurat <- readRDS("Seurat/VMH_IndependentAnalysis/VMHMerged_filtered3_Paperanalysis_30pcs.rds")
Male <- subset(x=mySeurat, subset=Hormone=="Intact")
Male <- subset(x=Male, idents=c("0"))
Primed <- subset(x=mySeurat, subset=Hormone=="Primed")
Primed <- subset(x=Primed, idents=c("0"))
Male.mat <- Male[["SCT"]]@scale.data
Male.mat <- as.matrix(Male.mat)
Primed.mat <- Primed[["SCT"]]@scale.data
Primed.mat <- as.matrix(Primed.mat)
signatures <- viper::rowTtest(x=Male.mat, y=Primed.mat)
#late_activity_matrix_Primed_vs_male <- viperSignature(eset=Primed.mat, ref=Male.mat, repos=TRUE, per=1000, cores=8)
#late_activity_matrix_scaled_Primed_vs_Male <- viper(late_activity_matrix_Primed_vs_male, regulon=regul, cores=8, verbose=F)
#head(late_activity_matrix_scaled_Primed_vs_Male)
#write.csv(late_activity_matrix_scaled_Primed_vs_Male, file="ViperOut/VMH/Vanillavipetest.csv")
signatures <- (qnorm(signatures$p.value/2, lower.tail=FALSE) * sign(signatures$statistic))[,1]
null_model <- viper::ttestNull(x=Male.mat, y=Primed.mat, per=1000, repos=TRUE, verbose=T, cores=8)
final_out <- viper::msviper(signatures, regul, null_model, verbose=T, cores=8)
saveRDS(final_out, file="ViperOut/VMH/Haniscripttest_1000permutecorrectcores.rds")
df_obj <- summary(final_out, length=(final_out$es$nes)) %>% as.data.frame() %>% arrange(FDR)
df_obj