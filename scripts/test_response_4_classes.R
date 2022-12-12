rm(list=ls())
require(randomForest)
require(dplyr)

test_beta <- readRDS(file = "../data/results/test_beta.Rds")
load(file="../out_100_probes_4_classes/rf.pred.RData")
beta <- t(test_beta)
probes <- rownames(importance(rf.pred, importance=1))
beta <- beta[, probes]
answer <- predict(rf.pred, beta, type = "prob")

response <- predict(rf.pred, beta, type="response")

test_1 <- colnames(answer)[apply(answer, 1, which.max)]
answer <- cbind(answer, test_1)
answer <- cbind(answer, levels(response)[response])
