#-----------------------------------------------------------------------------------
# training Random Forest classifier
#                                                                     
# Martin Sill
# m.sill@dkfz.de                                                                  
# 
# Note, in this example we reduce the number of features to 20k probes by sd filtering before applying the random forest for 
# feature selection. To perform feature selection as described in the paper remove line 34,
# this will increase the computation time significantly.
#
# 2018-03-14 UTC
#------------------------------------------------------------------------------------                   
#options(max.print = 1000)
#options(stringsAsFactors = FALSE)
#options(scipen = 999)
rm(list=ls())

library(randomForest)
library(doParallel)
setwd("..")

ntrees <- 10000  # 10000 in the paper, here 500 to speed up the example
cores <- detectCores()-1
seed <- 180314
p <- 10000   

message("loading preprocessed data ...",Sys.time())
#load(file.path("results","betas.ba.RData"))
if (!exists("betas"))
  betas <- as.data.frame(readRDS("./data/results/meth_combat_beta.Rds"))
if (!exists("meta_data"))
  meta_data <- read.table(file = "./data/meta_data/training_meta_data.txt", sep = "\t", header = T)

message("performing variable selection ...",Sys.time())
source(file.path("R","train.R"))
y <- as.factor(meta_data$CC_Epi_newLRO)

# sd pre filtering to 20k probes, to speed up the example
#betas <- betas[,order(-apply(betas,2,sd))]
betas <- t(betas)

set.seed(seed,kind ="L'Ecuyer-CMRG") 
message("seed: ",seed)
message("cores: ",cores)
message("ntrees: ",ntrees)  
message("n: ",nrow(betas))
message("p: ",ncol(betas))  

rf.varsel <- rfp(betas,
                 y,
                 mc=cores,
                 mtry = floor(sqrt(nrow(betas))),
                 ntree=ntrees,
                 tree_total=ntrees,
                 sampsize=rep(min(table(y)),length(table(y))),
                 importance=TRUE)

# get permutation variable importance
#imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1]
imp.meandecrease <- importance(rf.varsel, type=1)

# save selection forest
save(rf.varsel,file=file.path("results","varsel.RData"))
rm(rf.varsel)

# reduce data matrix
or <- order(imp.meandecrease,decreasing=T)
betasy <- betas[,or[1:p]]

gc()

message("finished ...",Sys.time())

message("training classifier ...",Sys.time())

message("single core")
message("ntrees: ",ntrees)  
message("n: ",nrow(betasy))
message("p: ",ncol(betasy))

rf.pred <- randomForest(betasy,
                        y,
                        #mc=cores,
                        ntree=ntrees,
                        #strata=y,
                        mtry=sqrt(ncol(betasy)),
                        sampsize=rep(min(table(y)),length(table(y))),
                        #proximity=TRUE,
                        #oob.prox=TRUE,
                        importance=TRUE,
                        #keep.inbag=TRUE,
                        do.trace=500,
                        seed=seed
)

message("finished ...",Sys.time())

save(rf.pred,file=file.path("results","rf.pred.RData"))

message("finished ...",Sys.time())
