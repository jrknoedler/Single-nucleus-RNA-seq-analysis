#!/usr/bin/env Rscript
data <- read.table("BNST_Compliliedforscript.txt", header=TRUE)
data <- data.frame(data)
type1 <- nrow(data[data$MvPfc > 0.5849 & data$MvPfdr < 0.05 & data$MvUfdr >0.05 & data$PvUfdr >0.05])

type2 <- nrow(data[data$MvPfdr > 0.05 & data$MvUfc > 0.5849 & data$MvUfdr <0.05 & data$PvUfdr >0.05])

type3 <- nrow(data[data$MvPfc > 0.5849 & data$MvPfdr <0.05 & data$MvUfc > 0.5849 & data$MvUfdr <0.05 & data$PvUfdr >0.05])

type4 <- nrow(data[data$MvPfc > 0.5849 & data$MvPfdr <0.05 & data$MvUfc > 0.5849 & data$MvUfdr <0.05 & data$PvUfc > 0.5849 & data$PvUfdr < 0.05])

type5 <- nrow(data[data$MvPfc > 0.5849 & data$MvPfdr <0.05 & data$MvUfc > 0.5849 & data$MvUfdr <0.05 & data$PvUfc < -0.5849 & data$PvUfdr < 0.05])

type6 <- nrow(data[data$MvPfc > 0.5849  & data$MvPfdr <0.05 & data$MvUfc > 0.5849 & data$MvUfdr <0.05 & data$PvUfc > 0.5849 & data$PvUfdr <0.05])

type7 <- nrow(data[data$MvPfc > 0.5849 & data$MvPfdr <0.05 & data$MvUfdr >0.05 & data$PvUfc < -0.5849 & data$PvUfdr <0.05])

type8 <- nrow(data[data$MvPfc < -0.5849 & data$MvPfdr <0.05 & data$MvUfdr >0.05 & data$PvUfdr >0.05])

type9 <- nrow(data[data$MvPfdr > 0.05 & data$MvUfdr >0.05 & data$PvUfc > 0.5849 & data$PvUfdr <0.05])

type10 <- nrow(data[data$MvPfc < -0.5849 & data$MvPfdr <0.05 & data$MvUfdr > 0.05 & data$PvUfc > 0.5849 & data$PvUfdr <0.05])

type11 <- nrow(data[data$MvPfc< -0.5849 & data$MvPfdr <0.05 & data$MvUfc > 0.5849 & data$MvUfdr < 0.05 & data$PvUfc > 0.5849 & data$PvUfdr < 0.05])

type12 <- nrow(data[data$MvPfc < -0.5849 & data$MvPfdr <0.05 & data$MvUfc < -0.5849 & data$MvUfdr < 0.05 & data$PvUfc > 0.5849 & data$PvUfdr < 0.05])

type13 <- nrow(data[data$MvPFC < -0.5849 & data$MvPfdr <0.05 & data$MvUfc < -0.5849 & data$MvUfdr < 0.05 & data$PvUfdr > 0.05])

type14 <- nrow(data[data$MvPfdr	> 0.05 & data$MvUfc < -0.5849 & data$MvUfdr <0.05 & data$PvUfdr > 0.05])

type15 <- nrow(data[data$MvPfdr	> 0.05 & data$MvUfdr > 0.05 & data$PvUfc < -0.5849 & data$PvUfdr < 0.05])

type16 <- nrow(data[data$MvPfc < -0.5849 & data$MvPfdr < 0.05 & data$MvUfdr > 0.05 & data$PvUfc < -0.5849 & data$PvUfdr < 0.05])

type17 <- nrow(data[data$MvPfc > 0.5849 & data$MvPfdr < 0.05 & data$MvUfc > 0.5849 & data$MvUfdr < 0.05 & data$PvUfc < -.05849 & data$PvUfdr < 0.05])

type18 <- nrow(data[data$MvPfc < -0.5849 & data$MvPfdr < 0.05 & data$MvUfc < -0.5849 & data$MvUfdr < 0.05 & data$PvUfc < -0.5849 & data$PvUfdr < 0.05])

compiled <- rbind(type1, type2, type3, type4, type5,type6,type7,type8,type9,type10,type11,type12,type13,type14,type15,type16,type17,type18)
write.csv(compiled, file="compiled_forceK.csv")