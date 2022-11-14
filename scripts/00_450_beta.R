rm(list = ls())
library(DMRcate)
library(dplyr)
library(ChAMP)
library(minfi)
library(readr)
library(base)
source("ChAMP_functions_Lionel_Philipp_220727.R")
#setwd("..")

# 450k pre-process
path_450k <- "./data/raw_data/training_450K.Rds"
meth_450.rgSet <- readRDS(path_450k)

## run NOOB 
meth_450 <- champ.load_extended(rgSet_object = meth_450.rgSet, 
                                sampleSheet = meth_450.sampleSheet,
                                method = "minfi", 
                                filterDetP = T, # include low p value probes
                                filterBeads = T, # include probes detected in few beads
                                beadCutoff = 1,
                                detPcut = 1,
                                arraytype = "450K", # set accorgingly
                                preproc = "Noob", 
                                dyeMethod = "single", 
                                dataToInclude = c("loadQC", "mset"),
                                force = F)
rm(meth_450.rgSet)

#remove legacy probes from 450k array
EPIC_legacy_probes <-  read_delim("./data/legacy_probes/MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2).txt.gz",
                                  delim = "\t",
                                  show_col_types = F)
meth_450$mset <- meth_450$mset[!(rownames(meth_450$mset) %in% EPIC_legacy_probes$TargetID), ]

## convert beta values for 450k sample
row_names_450 <- sort(rownames(meth_450$mset))
meth_450.beta <- getBeta(meth_450$mset, "Illumina")[row_names_450, ]
rm(meth_450, row_names_450)

mainDir <- ".."
subDir <- "pre_process"
dir.create(file.path(mainDir, subDir))
save(meth_450.beta,file=file.path("..","pre_process","00_450_beta.RData"))
