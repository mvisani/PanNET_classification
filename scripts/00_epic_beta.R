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

#load epic data
path_epic <- "./data/raw_data/training_EPIC.Rds"
meth_epic.rgSet <-  readRDS(path_epic)

meth_epic <- champ.load_extended(rgSet_object = meth_epic.rgSet, 
                                 sampleSheet = meth_epic.sampleSheet,
                                 method = "minfi", 
                                 filterDetP = T, # include low p value probes
                                 filterBeads = T, # include probes detected in few beads
                                 beadCutoff = 1,
                                 detPcut = 1,
                                 arraytype = "EPIC", # set accorgingly
                                 preproc = "Noob", 
                                 dyeMethod = "single", 
                                 dataToInclude = c("loadQC", "mset"),
                                 force = F)
rm(meth_epic.rgSet)

#Removing legacy probes. 
EPIC_legacy_probes <-  read_delim("./data/legacy_probes/MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2).txt.gz",
                                  delim = "\t",
                                  show_col_types = F)
meth_epic$mset <- meth_epic$mset[!(rownames(meth_epic$mset) %in% EPIC_legacy_probes$TargetID), ]
if(exists("EPIC.manifest.hg19"))
  rm(EPIC.manifest.hg19)


# convert to beta values (meth / unmeth + meth)
row_names_epic <- sort(rownames(meth_epic$mset))
meth_epic.beta <- getBeta(meth_epic$mset, "Illumina")[row_names_epic, ]


rm(meth_epic, row_names_epic) #remove data that we don't need anymore... (to test)
gc()

mainDir <- ".."
subDir <- "pre_process"
dir.create(file.path(mainDir, subDir))
save(meth_epic.beta,file=file.path("..","pre_process","00_epic_beta.RData"))

