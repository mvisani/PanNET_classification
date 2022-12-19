#rm(list=ls())
library(factoextra)
library(cluster)
require(readxl)

df <- readRDS(file="../data/results/meth_combat_beta.Rds")

alpha_vs_beta <-  read_xlsx(path="../data/raw_data/DiDomenico_2020_T6_alpha_like_v_beta_like.xlsx",
                            skip = 2)
alpha_vs_beta <- alpha_vs_beta[,1]
colnames(alpha_vs_beta) <- "probes"

alpha_vs_intermediate <- read_xlsx(path = "../data/raw_data/DiDomenico_2020_T7_alpha_like_v_intermediate.xlsx",
                                   skip = 2)
alpha_vs_intermediate <- alpha_vs_intermediate[,1]
colnames(alpha_vs_intermediate) <- "probes"

intermediate_vs_beta <- read_xlsx(path = "../data/raw_data/DiDomenico_2020_T8_intermediate_v_beta_like.xlsx",
                                  skip = 2)
intermediate_vs_beta <- intermediate_vs_beta[,1]
colnames(intermediate_vs_beta) <- "probes"

diff_meth_prob <- unique(rbind(intermediate_vs_beta, alpha_vs_intermediate, alpha_vs_beta))
diff_meth_prob <- intersect(diff_meth_prob$probes, rownames(df))
df <- df[diff_meth_prob,]
df <- t(df)

#gap_stat <- clusGap(df, FUN = hcut, K.max = 15, B = 20)
#fviz_gap_stat(gap_stat)

d <- dist(df, method = "euclidean")
final_clust <- hclust(d, method = "ward.D2" )
groups <- as.data.frame(cutree(final_clust, k=5))
table(groups)
save(groups, file = "../results/k_means_clustering.RData")