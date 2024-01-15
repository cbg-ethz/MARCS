#!/usr/bin/env Rscript

# Network-based clustering of mutational and covariate data
library(graphClust)

# load data
all_data <- read.csv("../data/undivided_binary_only_matrix_8k.csv")#[,c(1:65,71:74,77:78)]

# remove survival data, classifications and patient id 
index_remove <- c(which(colnames(all_data)=="ID"), which(colnames(all_data)=="OS"), 
                  which(colnames(all_data)=="OS_STATUS"), which(colnames(all_data)=="IPSSR_ELN"), 
                  which(colnames(all_data)=="WHO_2016"), which(colnames(all_data)=="WHO_2022"),
                  which(colnames(all_data)=="ICC"), which(colnames(all_data)=="ELN2022_IPSSM"))
mut_cov_data <- all_data[,-index_remove]
rownames(mut_cov_data) <- all_data$ID

# string_edgepmat <- as.matrix(read.table("../data/string-edgepmat.txt"))[-58,-58]
# blacklist_edgepmat <- as.matrix(read.table("../data/blacklist-edgepmat.txt"))[-58,-58]

# define covariates
n_cov <- 2

# clustering
set.seed(seednumber) # set seed

clusterResPlain <- get_clusters(mut_cov_data, k_clust = 12, n_bg = n_cov,
                                quick = FALSE, EMseeds = 1,
                                bdepar=list(chi = 0.5, edgepf = 8))
                                # edgepmat = string_edgepmat,
                                # blacklist = blacklist_edgepmat

saveRDS(clusterResPlain, paste0("euler_results/memberships_seed_", seednumber, ".rds"))


