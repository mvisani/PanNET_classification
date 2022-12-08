rm(list=ls())
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please specify your output directory", call.=FALSE)
}
outdir <- as.character(args[1])
dir.create(paste0("../", outdir),showWarnings = FALSE)

library(rmarkdown)
library(glmnet)
library(doParallel)
library(HandTill2001)


message("loading data ...",Sys.time())
if (!exists("betas"))
  betas <- readRDS("../data/results/meth_combat_beta.Rds")
if (!exists("meta_data"))
  meta_data <- read.table(file = "../data/meta_data/training_meta_data.txt", sep = "\t", header = T)

#load(file.path("results","Mset_filtered.RData"))
load(file.path("..",outdir,"nfolds.RData"))

for(i in 1:length(nfolds)){
  scores <- list() 
  idx <- list()
  for(j in 1:length(nfolds)){
    fname <- paste0("CVfold.",i,".",j,".RData")
    load(file.path("..",outdir,fname))
    scores[[j]] <- rf.scores
    idx[[j]] <- nfolds[[i]][[2]][[j]]$test
  }
  scores <- do.call(rbind,scores)
  idx <- unlist(idx)
  y <- meta_data$CC_Epi_newLRO[idx]         
  
  message("fitting calbriation model fold ",i," ...",Sys.time())
  # fit multinomial logistic ridge regression model
  suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                          alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))
  
  fname <- paste0("CVfold.",i,".",0,".RData")
  load(file.path("..", outdir,fname))
  
  message("calibrating raw scores fold ",i," ...",Sys.time())
  probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response"
                   ,s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
  
  
  err <- sum(colnames(probs)[apply(probs,1,which.max)] != meta_data$CC_Epi_newLRO[nfolds[[i]][[1]][[1]]$test])/length(nfolds[[i]][[1]][[1]]$test)
  
  message("misclassification error: ",err)
  
  fname_probs <- paste0("probsCVfold.",i,".",0,".RData")
  save(probs,file=file.path("..",outdir,fname_probs))
}

scores <- list()
idx <- list()
for(i in 1:length(nfolds)){
  fname <- paste0("CVfold.",i,".",0,".RData")
  load(file.path("..", outdir,fname))
  scores[[i]] <- rf.scores
  idx[[i]] <- nfolds[[i]][[1]][[1]]$test
}
scores <- do.call(rbind,scores)

probl <- list()
for(i in 1:length(nfolds)){
  fname <- paste0("probsCVfold.",i,".",0,".RData")
  load(file.path("..", outdir,fname))
  probl[[i]] <- probs
}
probs <- do.call(rbind,probl)


idx <- unlist(idx)
y <- meta_data$CC_Epi_newLRO[idx] 

ys <- colnames(scores)[apply(scores,1,which.max)]
yp <- colnames(probs)[apply(probs,1,which.max)]

errs <- sum(y!=ys)/length(y)
errp <- sum(y!=yp)/length(y)

message("overall misclassification error scores: ",errs)
message("overall misclassification error calibrated: ",errp)

message("fitting final calibration model ...",Sys.time())

suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                        alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))

save(cv.calfit,file=file.path("..", outdir,"calfit.RData"))

save(scores,probs,y,ys,yp,file=file.path("..", outdir,"CVresults.RData"))

#message("generating report ...",Sys.time())
#rmarkdown::render("CVresults.Rmd")
message("finished ...",Sys.time())