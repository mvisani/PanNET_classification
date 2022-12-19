rm(list = ls())
require(ChAMP)
library(DMRcate)
library(dplyr)
library(ChAMP)
library(minfi)
library(readr)
library(base)
library(doParallel)
source("ChAMP_functions_Lionel_Philipp_220727.R")

df_epic <- readRDS(file="../data/raw_data/training_EPIC.Rds")
p_val_epic <- as.data.frame(detectionP(df_epic))
rm(df_epic)
names_failed_epic <- as.data.frame(rownames(p_val_epic)
                                   [apply(p_val_epic, 1, function(x) sum(x>= 0.01))>=ceiling(0.04*dim(p_val_epic)[2])])
colnames(names_failed_epic) <- "probes"

df_450 <- readRDS(file="../data/raw_data/training_450K.Rds")
p_val_450 <- as.data.frame(detectionP(df_450))
rm(df_450)
names_failed_450 <- as.data.frame(rownames(p_val_450)
                                  [apply(p_val_450, 1, function(x) sum(x>= 0.01))>=ceiling(0.04*dim(p_val_450)[2])])
colnames(names_failed_450) <- "probes"

names_failed <- as.data.frame(unique(rbind(names_failed_epic, names_failed_450)))

require(randomForest)
load(file="../out_100_probes/rf.pred.RData")
probes <- as.data.frame(rownames(importance(rf.pred, type=1)))
colnames(probes) <- "probes"
length(intersect(probes$probes, names_failed$probes))/length(probes$probes)

