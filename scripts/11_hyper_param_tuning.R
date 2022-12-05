#rm(list=ls())

library(randomForest)
library(minfi)
library(limma)
#setwd("..")

ntrees <- c(100, 500, 1000, 5000, 10000)
cores <- detectCores() - 1
seed <- 180314
p <- c(10, 100, 1000, 10000)
folds <- 3

if (!exists("betas"))
  betas <- readRDS("../data/results/meth_combat_beta.Rds")
if (!exists("meta_data"))
  meta_data <- read.table(file = "../data/meta_data/training_meta_data.txt",
                          sep = "\t", header = T)

y <- as.factor(meta_data$CC_Epi_newLRO)
#betas <- t(betas)

load(file = "../results/varsel.RData")
imp.meandecrease <- importance(rf.varsel, type=1)
or <- order(imp.meandecrease,decreasing=T)

par(mfrow=c(4, 5),
    mar = c(2, 2, 2, 2))
for (probes in p) {
  for (j in ntrees) {
    cv <- rfcv(betas[,or[1:probes]],
         y,
         cv.fold = folds, 
         ntree = j,
         seed = seed)
    with(cv, plot(n.var,
                  error.cv,
                  type="o",
                  lwd=2,
                  log="x",
                  ylim=c(0,1),
                  xlab="n.probes",
                  ylab="error.cv", 
                  main= paste0("Top ", probes, " probes and delopped with ", j, " trees."),
                  cex.main = 0.7))
    abline(h=0.2, lty=2)
  }
}
