rm(list=ls())
require(randomForest)
require(dplyr)

test_beta <- readRDS(file = "../data/results/test_beta.Rds")
meta_data <- read.table(file="../data/meta_data/test_EPIC_meta_data.txt", sep = "\t", header = T)
load(file="../out_100_probes/rf.pred.RData")
beta <- t(test_beta)
probes <- rownames(importance(rf.pred, importance=1))
beta <- beta[, probes]
answer <- predict(rf.pred, beta, type = "prob")

prediction <- colnames(answer)[apply(answer, 1, which.max)]
answer <- cbind(answer, prediction)


CC_epi <- meta_data$CC_Epi_newG1G2_plus_met
answer <- cbind(answer, CC_epi)

answer <- as.data.frame(answer)
answer["Sample_Name"] <- rownames(answer)

require(dplyr)
doubled <- meta_data %>% group_by(aPmP_number) %>% filter(n()>1)

doubled <- answer[doubled$Sample_Name, ]
doubled$prediction==doubled$CC_epi
