#this scirpt joins the datasets coming from the 2 different technologies for the training data.
# Normally it shouldn't be used anymore. 
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

load("../pre_process/00_epic_beta.RData")
load("../pre_process/00_450_beta.RData")

# Combine 450K and EPIC beta matrices (intersection)
row_names <- intersect(rownames(meth_epic.beta), rownames(meth_450.beta))
meth.beta <- cbind(meth_450.beta[row_names, ], meth_epic.beta[row_names,])

n_cores <- detectCores() - 1
#Normalize dataset
meth.norm <- champ.norm(meth.beta, 
                        arraytype = "450K", 
                        cores = n_cores,
                        resultsDir = "../data/results")
saveRDS(meth.norm, file="../data/results/meth_normalized.Rds", compress = T)
if(exists("meth.beta"))
  rm(meth.beta)