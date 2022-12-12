rm(list=ls())
library(cluster)
library(factoextra)
require(fclust)
require(e1071)
require(corrplot)

betas <- readRDS(file = "../data/results/meth_combat_beta.Rds")
meta_data <- read.table(file = "../data/meta_data/training_meta_data_new.txt", sep = "\t", header = T)
y <- as.factor(meta_data$CC_Epi_newLRO)
probvar <- apply(betas, 1, var)
probvar <- order(probvar, decreasing = T)
betas <- betas[probvar[1:100000], ]
betas <- t(betas)

k <- kmeans(betas, centers = 6, iter.max = 100)
fviz_cluster(k , data = betas)

cm_3 <- cmeans(betas, centers = 3, iter.max = 2000, m=3)
cm_4 <- cmeans(betas, 4, iter.max = 200)
cm_5 <- cmeans(betas, 5, iter.max = 200)
cm_6 <- cmeans(betas, 6, iter.max = 200)

fviz_cluster(list(data = betas, cluster=cm_3$cluster), 
             palette = "jco",
             ggtheme = theme_minimal())
corrplot(cm_3$membership, tl.cex = 0.2)
