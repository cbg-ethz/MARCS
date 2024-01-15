# get the best performing seed and corresponding DAGs + memberships
library(clustNet)
library(graphClust)
rm(list=ls())

assignprogressList <- list()
# seeds <- which(sapply(1:100, function(x) length(readRDS(paste0("euler_results/memberships_seed_",x,".rds"))$DAGs))==9)
seeds <- c(1:100)
for (i in seeds){
  assignprogressList[[i]] <- readRDS(paste0("euler_results/memberships_seed_",i,".rds"))$assignprogress
}

myseed <- function (assignprogress) 
{
  likelihoodvector <- unlist(lapply(assignprogress, function(x) x[length(x)]))
  EMseeds <- 1:length(assignprogress)
  EMn <- length(EMseeds)
  so <- order(likelihoodvector, decreasing = TRUE)
  bestseed <- EMseeds[so][1]
  worstseed <- EMseeds[so][EMn]
  return(bestseed)
}

#best_seed <- clustNet:::getBestSeed(assignprogressList)
best_seed <- myseed(assignprogressList)
best_res <- readRDS(paste0("euler_results/memberships_seed_",seeds[best_seed],".rds"))

saveRDS(best_res,"./AMLonly_6clus.rds")

