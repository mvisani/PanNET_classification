#!/usr/bin/env Rscript
rm(list=ls())
library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="raw data file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name", metavar="character"),
  make_option(c("--outdir"), type = "character", default = "../data/results", metavar = "character",
              help="output directory"),
  make_option(c("--legacy_probes"), type = "character", metavar = "character",
              default = "../data/legacy_probes/MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2).txt.gz",
              help = "remove the legacy probes of specified file [default= %default]"),
  make_option(c("-m", "--meta_data"), type = "character", metavar = "character", 
              default = NULL,
              help="meta data file")
); 

opt_parser <- OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (is.null(opt$meta_data)){
  print_help(opt_parser)
  stop("At least one argument must be supplied for the meta data (input file).n", call.=FALSE)
}

library(DMRcate)
library(dplyr)
library(ChAMP)
library(minfi)
library(readr)
library(base)
library(doParallel)
source("ChAMP_functions_Lionel_Philipp_220727.R")

path_epic <- opt$file
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
EPIC_legacy_probes <-  read_delim(opt$legacy_probes,
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

dir.create(opt$outdir)
n_cores <- detectCores() - 1
#Normalize dataset
meth.norm <- champ.norm(meth_epic.beta, 
                        arraytype = "450K", 
                        cores = n_cores,
                        resultsDir = opt$outdir)


message("removing batch effect with ComBat")
if(!exists("meta_data"))
  meta_data <- read.table(opt$meta_data, sep = "\t", header = T)

meth.norm <- meth.norm[, meta_data$Sample_Name]

meth_cb_model <- model.matrix(~ 1 , data = meta_data)

require(sva)
meth_combat <- ComBat(ENmix::B2M(meth.norm),
                      batch = meta_data$Slide,
                      mod = meth_cb_model,
                      BPPARAM = bpparam("SerialParam"))

saveRDS(meth_combat, file = paste0(opt$outdir, opt$out, "_M.Rds"), compress = T)

meth_combat_beta <- ENmix::M2B(meth_combat)
saveRDS(meth_combat_beta, file=paste0(opt$outdir, opt$out,"_B.Rds"))
message("Done. Have a great day :) ")