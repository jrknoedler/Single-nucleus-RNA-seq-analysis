#!/usr/bin/env Rscript

sampleID <- "MalePOA"
directory <- "MalePOAIntact/outs/filtered_feature_bc_matrix"
output <- "Seurat/POA_IndependentAnalysis/POA_Viper_"

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
regul
#write.csv(regul, file="regulon.csv")
mySeurat <- readRDS("Seurat/POA_IndependentAnalysis/POA_Fullmerge_Kiss1opt_Filtered3_40pcs_res1.2.rds")
mySeurat <- SCTransform(mySeurat, vars.to.regres="orig.ident", return.only.var.genes=FALSE)

for (i in 0:34){
try({
ident <- i

Male <- subset(x=mySeurat, subset=Hormone=="Intact")
Esr1 <- subset(x=mySeurat, idents=c("10","16","20","26","7","11","3","21","14","21","4","19","28"))
Primed <- subset(x=mySeurat, subset=Hormone=="Unprimed")
Ref <- subset(x=mySeurat, idents=c("10","16","20","26","7","11","3","21","14","21","4","19","28"), invert=TRUE)
Ref.mat <- Ref[["SCT"]]@scale.data
Ref.mat <- as.matrix(Ref.mat)
dim(Ref.mat)
Esr1.mat <- Esr1[["SCT"]]@scale.data
Esr1.mat <- as.matrix(Esr1.mat)
signatures <- viper::rowTtest(x=Esr1.mat, y=Ref.mat)
#late_activity_matrix_Primed_vs_male <- viperSignature(eset=Primed.mat, ref=Male.mat, repos=TRUE, per=1000, cores=8)
#late_activity_matrix_scaled_Primed_vs_Male <- viper(late_activity_matrix_Primed_vs_male, regulon=regul, cores=8, verbose=F)
#head(late_activity_matrix_scaled_Primed_vs_Male)
#write.csv(late_activity_matrix_scaled_Primed_vs_Male, file="ViperOut/VMH/Vanillavipetest.csv")
signatures <- (qnorm(signatures$p.value/2, lower.tail=FALSE) * sign(signatures$statistic))[,1]
null_model <- viper::ttestNull(x=Esr1.mat, y=Ref.mat, per=1000, repos=TRUE, verbose=T)
final_out <- viper::msviper(signatures, regul, null_model, verbose=T)
saveRDS(final_out, file="ViperOut/VMH/Haniscripttest_1000permutecorrectcores.rds")
df_obj <- summary(final_out, length=(final_out$es$nes)) %>% as.data.frame() %>% arrange(FDR)
df_obj
pdf(file=paste0(output,"mrsplot1Esr1.pdf"))
plot(final_out, cex=0.7)
dev.off()
mrs <- ledge(final_out)
summary(mrs)