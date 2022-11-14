rm(list = ls())
library(DMRcate)
library(dplyr)
library(ChAMP)
library(minfi)
library(readr)
library(base)
library(doParallel)
source("ChAMP_functions_Lionel_Philipp_220727.R")
#setwd("..")

if(!exists("meth.norm"))
  meth.norm <- readRDS("../data/results/meth_normalized.Rds")
if(!exists("meta_data"))
  meta_data <- read.table("../data/meta_data/training_meta_data.txt",sep = "\t", header = T)

meth.norm <- meth.norm[, meta_data$Sample_Name]

meth_cb_model <- model.matrix(~ 1 , data = meta_data)

require(sva)
meth_combat <- ComBat(ENmix::B2M(meth.norm),
                      batch = meta_data$Slide,
                      mod = meth_cb_model,
                      BPPARAM = bpparam("SerialParam"))

saveRDS(meth_combat, file = "../data/results/meth_combat_M.Rds", compress = T)

meth_combat_beta <- ENmix::M2B(meth_combat)
saveRDS(meth_combat_beta, file="../data/results/meth_combat_beta.Rds")