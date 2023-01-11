# output an image showing the model fine tuned. One should specifiy the number of probes and trees
# and also if one wants to check for 4 classes model or 5 classes model. 
rm(list=ls())
library(randomForest)
library(minfi)
library(limma)
library(doParallel)

cores <- detectCores() - 1
#ntrees <- c(100, 500, 1000, 5000, 10000)
ntrees <- c(500)
seed <- 180314
#p <- c(10, 100, 1000, 10000)
p <- c(10, 25, 50, 100, 200, 500)
folds <- 3

if (!exists("betas"))
  betas <- readRDS("../data/results/meth_combat_beta.Rds")
if (!exists("meta_data"))
  meta_data <- read.table(file = "../data/meta_data/training_meta_data_new.txt",
                          sep = "\t", header = T)
y <- as.factor(meta_data$Four_classes)
source(file.path("R","makefolds.R"))
source(file.path("R", "train.R"))
nfolds <- makenestedfolds(y,folds)

betas <- t(betas)
rf.varsel <- rfp(betas,
                 y,
                 mc=cores,
                 mtry = floor(sqrt(nrow(betas))),
                 ntree=1000,
                 sampsize=rep(min(table(y)),length(table(y))),
                 importance=TRUE, 
                 seed = seed)


imp.meandecrease <- importance(rf.varsel, type=1)
or <- order(imp.meandecrease,decreasing=T)

png(filename = "../results/n_probe_selection_4_classes.png", width = 2000, height = 1538)
#par(mfrow=c(length(p), length(ntrees)), mar = c(2, 2, 2, 2))
par(mfrow=c(2,3))
for (probes in p) {
  for (j in ntrees) {
    plot(NA, main= paste0("Top ", probes, " probes with ", j, " trees."),
         xlab="n.probes",
         ylab="error.cv",
         ylim=c(0,1),
         xlim=c(1, probes),
         cex.main = 1.5,
         log="x"
         )
    
    for (i in 1:length(nfolds)){
      samp <- nfolds[[i]][[1]][[1]][["train"]]
      cv <- rfcv(betas[samp,or[1:probes]],
                 y[samp],
                 cv.fold = folds, 
                 ntree = j,
                 do.trace=T,
                 seed = seed)
      with(cv, lines(n.var,
                    error.cv,
                    type="o",
                    lwd=2,
                    col=i))
    }
    abline(h=0.2, lty=2)
  }
}
dev.off()