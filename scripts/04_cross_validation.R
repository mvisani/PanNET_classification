# Perform cross validation on the training data. 
# The script will perform inner fold cross validations on model. 
# for further explanation see : http://www.nature.com/doifinder/10.1038/nature26000 
rm(list=ls())
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please specify your output directory as first argument and the number of probes as second argument", call.=FALSE)
}

library(randomForest)
library(doParallel)
library(minfi)
library(limma)

ntrees <- 500 
cores <- detectCores()
seed <- 180314
folds <- 3

outdir <- as.character(args[1])
dir.create(paste0("../", outdir),showWarnings = FALSE)
p <- as.double(args[2])

if (!exists("betas"))
  betas <- readRDS("../data/results/meth_combat_beta.Rds")
if (!exists("meta_data"))
  meta_data <- read.table(file = "../data/meta_data/training_meta_data.txt", sep = "\t", header = T)

y <- as.factor(meta_data$CC_Epi_newLRO)
batch <- as.factor(meta_data$Technology)

source(file.path("R","makefolds.R"))
source(file.path("R","train.R"))
source(file.path("R","calculateCVfold.R"))
source(file.path("R","batchadjust.R"))

if(!file.exists(file.path("..", outdir,"nfolds.RData"))){
  nfolds <- makenestedfolds(y,folds)
  save(nfolds,file=file.path("..",outdir,"nfolds.RData"))
}
load(file.path(".." ,outdir,"nfolds.RData"))

message("performing nested CV ...", Sys.time())
message("check minimal class sizes for inner training loops")

# check minimal class sizes for inner training loops
minclasssize <- matrix(0,ncol=length(nfolds),nrow=length(nfolds))
for(i in 1:length(nfolds)){
  for(j in 1:length(nfolds))
    minclasssize[i,j]  <- min(table(y[nfolds[[i]][[2]][[j]]$train]))
}
colnames(minclasssize) <- paste0("innfold",1:folds)
rownames(minclasssize) <- paste0("fold",1:folds)
print(minclasssize)

for(K in 1:folds){
  
  for(k in 0:folds){
    
    if(k>0){  message("calculateing fold ",K,".",k,"  ...",Sys.time())
      fold <- nfolds[[K]][[2]][[k]]
    }else{
      message("calculateing outer fold ",K,"  ...",Sys.time())
      fold <- nfolds[[K]][[1]][[1]]
    }
    
    rf.scores <- calcultateCVfold(betas,y,batch,fold,p,cores,ntrees)
    
    fname <- paste("CVfold",K,k,"RData",sep=".")
    save(rf.scores,file=file.path("..", outdir,fname))
    
    rm(rf.scores)
    gc()
  }
}
message("finished ...",Sys.time())
