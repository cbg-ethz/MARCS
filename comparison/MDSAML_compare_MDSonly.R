library(graphClust)
library(scales)
library(ggplot2)
library(survival)
library(dplyr)
# Load required libraries
library(survival)
library(survminer)

# # Simulate binary data from 3 clusters
# k_clust <- 9
# all_data <- read.csv("../data/undivided_binary_matrix.csv")
# 
# clustered_data <- all_data[,c(3:66,69:74)]
# 
# # put age and sex in last row
# clustered_data[,c(1:64,67:70,65:66)]
# clustered_data$AGE <- clustered_data$AGE-1
# 
# #store results
# cluster_res <- readRDS("../results/custer_res.rds")

# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)

# Load classification by mutation profile (cluster assignment) # CLUSTERS 

#cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")
cluster_res <- readRDS("./cluster_MDS/MDSonly_6clus.rds")
#clusterMDS <- readRDS("./MDSonly_6clus.rds")
#readRDS("../results/euler_memberships_8k_9clusters.rds")

# remove this line in the future (after package update on euler)
# cluster_res$clustermembership <- cluster_res$clustermembership[[1]]
# cluster_res <- readRDS("../results/custer_res.rds")

table(cluster_res$clustermembership)



# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
mut_cov_data <- read.csv("../data/undivided_binary_only_matrix_8k.csv")
# survdata <- survdata[survdata$ID %in% rownames(cluster_res$data),]
# # clumatrix <- data.frame(ID=rownames(cluster_res$data), group=cluster_res$clustermembership)
# # survdata <- merge(survdata,clumatrix, by = "ID")

# survdata$last_fu <- as.Date(survdata$last_fu, format = "%d.%m.%Y")
# survdata$firstdiagnosis_date <- as.Date(survdata$firstdiagnosis_date, format = "%d.%m.%Y")
# survdata$relapse_date <- as.Date(survdata$relapse_date, format = "%d.%m.%Y")
# survdata$type <- as.factor(survdata$)
survdata <- survdata[grep("MDS/AML",survdata$ICC),]

mut_cov_data <- mut_cov_data[mut_cov_data$ID %in% survdata$ID,]

# ids <- mutation_covariate_data$UPI[mutation_covariate_data$Dx %in% c("AML","MDS")]
# validata <- validata[validata$UPI %in% ids,]

# totdata <- merge(totdata, riskstrat[,c(1,85)], by="UPI")
# totdata <- totdata[totdata$UPI %in% ids,]

# define the function
get_classification <- function(cluster_results, data_classify){
  
  myData <- cluster_results$data
  k_clust <- length(cluster_results$DAGs)
  
  if(is.vector(data_classify)){
    data_classify <- t(as.data.frame(data_classify)) # when this is a single col entry
  }else{
    data_classify <- as.data.frame(data_classify)
  }
  
  # input is clustercenters
  clustercenters <- cluster_results$DAGs
  newallrelativeprobabs <- cluster_results$probs
  
  
  ## detect and adjust for missing data
  
  # Create two example matrices
  matrix_one <- myData
  
  # Create two example matrices
  matrix_two <- data_classify
  
  # Find shared variables
  shared_vars <- intersect(colnames(matrix_one), colnames(matrix_two))
  # cat("Shared variables:", shared_vars, "\n")
  
  # Find missing variables in matrix_one
  missing_vars_one <- setdiff(colnames(matrix_one), colnames(matrix_two))
  if(length(missing_vars_one)> 0){
    cat("Missing variables in cluster data:", missing_vars_one, "\n")
  }
  
  # Find missing variables in matrix_two
  missing_vars_two <- setdiff(colnames(matrix_two), colnames(matrix_one))
  if(length(missing_vars_two)> 0){
    cat("Missing variables in classification data:", missing_vars_two, "\n")
  }
  
  # Get column indices of shared variables in matrix_one
  if (length(shared_vars) > 0) {
    shared_columns_indices_one <- match(shared_vars, colnames(matrix_one))
    shared_columns_indices_two <- match(shared_vars, colnames(matrix_two))
    # adapt data to shared vars
    myData <- matrix_one[,shared_columns_indices_one]
    data_classify <- matrix_two[,shared_columns_indices_two]
    # adapt DAGs to shared vars
    for (ii in 1:k_clust){
      clustercenters[[ii]] <- clustercenters[[ii]][shared_columns_indices_one,shared_columns_indices_one]
    }
  } else {
    cat("No shared variables\n")
  }
  
  ## classification
  
  data_classify <- unname(data_classify)
  
  # number of background variables
  n_bg <- 2
  # total number of variables
  ss <- dim(myData)[1]
  nn <- dim(myData)[2]
  ss2 <- NROW(data_classify)
  #number of variables (without covariates)
  n<-nn-n_bg
  
  scoresagainstclusters<-matrix(ncol=k_clust,nrow=ss2)
  
  bdepar <- list(chi = 0.5, edgepf = 8)
  
  
  if (!all(myData < 2)){
    # cetegorical version
    score_type <- "bdecat"
  }else{
    # binary version
    score_type <- "bde"
  }
  
  allrelativeprobabs<-newallrelativeprobabs
  coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
  tauvec<-coltots/sum(coltots)
  
  parRes <- parallel::mclapply(1:k_clust, function(k) {
    if (n_bg>0){
      scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData),
                                         weightvector=allrelativeprobabs[,k],
                                         bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
    }else{
      scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData), edgepmat = edgepmat,
                                         weightvector=allrelativeprobabs[,k],
                                         bdepar=bdepar)
    }
    
    scorepar$n <- n # to avoid to scoring over background nodes
    # scoresagainstclusters[,k] <- BiDAG::scoreagainstDAGscoreagainstDAG(scorepar,clustercenters[[k]])
    
    # if (score_type=="bdecat"){
    #   scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], bdecatCvec = apply(myData, 2, function(x) length(unique(x))))
    # }else{
    #   scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]])
    # }
    
    if (score_type=="bdecat"){
      scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], bdecatCvec = apply(myData, 2, function(x) length(unique(x))), datatoscore = data_classify)
    }else{
      scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], datatoscore = data_classify)
    }
    
    
    scorepar$n <- n+n_bg # recet after scoring
    
    return(scoresagainstclusters[,k])
    
    # }, mc.cores = k_clust)
  })
  
  for (kk in 1:k_clust){
    scoresagainstclusters[,kk] <- parRes[[kk]]
  }
  
  # assign cluster
  newallrelativeprobabsnotau <- graphClust:::allrelativeprobs(scoresagainstclusters)
  newallrelativeprobabs <- graphClust:::relativeprobswithtau(newallrelativeprobabsnotau,tauvec)
  
  newclustermembership<-graphClust:::reassignsamples(newallrelativeprobabs)
  
  return(list("clustermembership"=newclustermembership,"allrelativeprobabs"=newallrelativeprobabs, "shared_vars"=shared_vars, "classified_data"=myData))
}

reclass <- get_classification(cluster_res, mut_cov_data)
survdata$group <- reclass$clustermembership
survdata$Cluster <- survdata$group
survdata$Cluster <- as.factor(survdata$Cluster)

p <- c()

for(i in unique(survdata$ICC)){
  subdata <- survdata[survdata$ICC== i,]
  group_count <- subdata %>% 
    group_by(Cluster) %>% 
    dplyr::summarise(n = n())
  
  filtered_groups <- group_count %>% 
    filter(n >= 10) %>% 
    pull(Cluster)
  
  # Create a new DataFrame containing only the rows where group appears 5 or more times
  subdata <- subdata %>% 
    filter(Cluster %in% filtered_groups)
  
  # merge clinical information and cluster membership
  clinical <- list()
  clinical$event <- c(subdata$OS_STATUS,subdata$OS_STATUS)
  clinical$time <- c(subdata$OS,subdata$OS)
  clinical$group <- c(subdata$ICC,subdata$Cluster)
  clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)
  
  # THESE NOW MATCH
  length(clinical$group)
  length(clinical$time)
  
  table(clinical$group)
  
  # Kaplan-Meier curve for groups
  
  nb.cols <- length(unique(cluster_res$clustermembership))
  
  # 
  # 
  library(viridis)
  #colours_clusters <- c("blue", "#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")
  colours_clusters <- c("blue", viridis(6))
  #plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
  #legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')
  mypal=colours_clusters[c(1, c(seq(1:9)+1)[seq(1,5) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])]
  dimensione=length(mypal)-1
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
  library("survminer")
  file<-paste0("../comparison/Supp_",gsub("MDS/AML ", "MDSAML_",i),"MDS_only.png")
  png(filename = file, width = 2300, height = 2000, res = 280)
  p[i] <- ggsurvplot(
    os,                     # survfit object with calculated statistics.
    data = clinical,             # data used to fit survival curves.
    #palette = colours_clusters,#colours_clusters, # personalized colours
    palette = alpha(mypal,c(1,rep(0.55, dimensione))), 
    # legend.labs = sort(as.character(unique(survdata$Cluster))), 
    risk.table = F,       # show risk table.
    pval = T,             # show p-value of log-rank test.
    conf.int = F,         # show confidence intervals for 
    # point estimates of survival curves.
    xlim = c(0,10),         # was c(0,14.4) 
    # survival estimates.
    xlab = "Time in years",   # customize X axis label.
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.fontsize = 3,
    pval.size =4,
    pval.coord = c(0, 0.15),
    title= paste0("OS - MDS only -", gsub(pattern = "_", replacement = " ", i), " (n=",length(clinical$event)/2,")"),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
    
  )
  print(p[i])
  
  dev.off()
}

pdf(file = "../analysis/figures/Supplementary_MDSAML_OS_MDSonly.pdf", width = 16, height = 10)
ggarrange(p$`MDS/AML MR-Cyto`,p$`MDS/AML MR-Gen`, p$`MDS/AML TP53`, p$`MDS/AML NOS`, ncol = 2, nrow = 2 , labels = letters[1:4])
dev.off() 



#### PFS

library(graphClust)
library(scales)
library(ggplot2)
library(survival)
library(dplyr)
# Load required libraries
library(survival)
library(survminer)


# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)

# Load classification by mutation profile (cluster assignment) # CLUSTERS 
#cluster_res <- readRDS("../results/euler_memberships_final_9.rds")
# cluster_res <- readRDS("./results/euler_memberships_final_26.rds")
cluster_res <- readRDS("./MDSonly_6clus.rds")
#clusterMDS <- readRDS("./MDSonly_6clus.rds")
#readRDS("../results/euler_memberships_8k_9clusters.rds")

# remove this line in the future (after package update on euler)
# cluster_res$clustermembership <- cluster_res$clustermembership[[1]]
# cluster_res <- readRDS("../results/custer_res.rds")

table(cluster_res$clustermembership)



# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
mut_cov_data <- read.csv("../data/undivided_binary_only_matrix_8k.csv")
# survdata <- survdata[survdata$ID %in% rownames(cluster_res$data),]
# # clumatrix <- data.frame(ID=rownames(cluster_res$data), group=cluster_res$clustermembership)
# # survdata <- merge(survdata,clumatrix, by = "ID")

# survdata$last_fu <- as.Date(survdata$last_fu, format = "%d.%m.%Y")
# survdata$firstdiagnosis_date <- as.Date(survdata$firstdiagnosis_date, format = "%d.%m.%Y")
# survdata$relapse_date <- as.Date(survdata$relapse_date, format = "%d.%m.%Y")
# survdata$type <- as.factor(survdata$)
# survdata$last_fu <- as.Date(survdata$last_fu, format = "%d.%m.%Y")
# survdata$firstdiagnosis_date <- as.Date(survdata$firstdiagnosis_date, format = "%d.%m.%Y")
# survdata$relapse_date <- as.Date(survdata$relapse_date, format = "%d.%m.%Y")
# survdata$type <- as.factor(survdata$)
survdata$WHO_2016 <- as.factor(survdata$WHO_2016)
survdata$WHO_2022 <- as.factor(survdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(survdata$IPSSR_ELN)
survdata$ELN2022_IPSSM <- as.factor(survdata$ELN2022_IPSSM)
pfsdata <- read.delim("../data/IPSSM_df_clinical.tsv")
pfsdata <- pfsdata[, c("ID","AMLt_YEARS","AMLt_STATUS","LFS_YEARS","LFS_STATUS")]
survdata <- merge(survdata,pfsdata, by = "ID")
survdata$WHO_2016 <- as.factor(survdata$WHO_2016)
survdata$WHO_2022 <- as.factor(survdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(survdata$IPSSR_ELN)
survdata$ICC <- as.factor(survdata$ICC)

survdata <- survdata[grep("MDS/AML",survdata$ICC),]

mut_cov_data <- mut_cov_data[mut_cov_data$ID %in% survdata$ID,]

# define the function
get_classification <- function(cluster_results, data_classify){
  
  myData <- cluster_results$data
  k_clust <- length(cluster_results$DAGs)
  
  if(is.vector(data_classify)){
    data_classify <- t(as.data.frame(data_classify)) # when this is a single col entry
  }else{
    data_classify <- as.data.frame(data_classify)
  }
  
  # input is clustercenters
  clustercenters <- cluster_results$DAGs
  newallrelativeprobabs <- cluster_results$probs
  
  
  ## detect and adjust for missing data
  
  # Create two example matrices
  matrix_one <- myData
  
  # Create two example matrices
  matrix_two <- data_classify
  
  # Find shared variables
  shared_vars <- intersect(colnames(matrix_one), colnames(matrix_two))
  # cat("Shared variables:", shared_vars, "\n")
  
  # Find missing variables in matrix_one
  missing_vars_one <- setdiff(colnames(matrix_one), colnames(matrix_two))
  if(length(missing_vars_one)> 0){
    cat("Missing variables in cluster data:", missing_vars_one, "\n")
  }
  
  # Find missing variables in matrix_two
  missing_vars_two <- setdiff(colnames(matrix_two), colnames(matrix_one))
  if(length(missing_vars_two)> 0){
    cat("Missing variables in classification data:", missing_vars_two, "\n")
  }
  
  # Get column indices of shared variables in matrix_one
  if (length(shared_vars) > 0) {
    shared_columns_indices_one <- match(shared_vars, colnames(matrix_one))
    shared_columns_indices_two <- match(shared_vars, colnames(matrix_two))
    # adapt data to shared vars
    myData <- matrix_one[,shared_columns_indices_one]
    data_classify <- matrix_two[,shared_columns_indices_two]
    # adapt DAGs to shared vars
    for (ii in 1:k_clust){
      clustercenters[[ii]] <- clustercenters[[ii]][shared_columns_indices_one,shared_columns_indices_one]
    }
  } else {
    cat("No shared variables\n")
  }
  
  ## classification
  
  data_classify <- unname(data_classify)
  
  # number of background variables
  n_bg <- 2
  # total number of variables
  ss <- dim(myData)[1]
  nn <- dim(myData)[2]
  ss2 <- NROW(data_classify)
  #number of variables (without covariates)
  n<-nn-n_bg
  
  scoresagainstclusters<-matrix(ncol=k_clust,nrow=ss2)
  
  bdepar <- list(chi = 0.5, edgepf = 8)
  
  
  if (!all(myData < 2)){
    # cetegorical version
    score_type <- "bdecat"
  }else{
    # binary version
    score_type <- "bde"
  }
  
  allrelativeprobabs<-newallrelativeprobabs
  coltots<-colSums(allrelativeprobabs) + bdepar$chi # add prior to clustersizes
  tauvec<-coltots/sum(coltots)
  
  parRes <- parallel::mclapply(1:k_clust, function(k) {
    if (n_bg>0){
      scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData),
                                         weightvector=allrelativeprobabs[,k],
                                         bdepar=bdepar, bgnodes=(n+1):(n+n_bg))
    }else{
      scorepar <- BiDAG::scoreparameters(score_type,as.data.frame(myData), edgepmat = edgepmat,
                                         weightvector=allrelativeprobabs[,k],
                                         bdepar=bdepar)
    }
    
    scorepar$n <- n # to avoid to scoring over background nodes
    # scoresagainstclusters[,k] <- BiDAG::scoreagainstDAGscoreagainstDAG(scorepar,clustercenters[[k]])
    
    # if (score_type=="bdecat"){
    #   scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], bdecatCvec = apply(myData, 2, function(x) length(unique(x))))
    # }else{
    #   scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]])
    # }
    
    if (score_type=="bdecat"){
      scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], bdecatCvec = apply(myData, 2, function(x) length(unique(x))), datatoscore = data_classify)
    }else{
      scoresagainstclusters[,k] <- BiDAG::scoreagainstDAG(scorepar,clustercenters[[k]], datatoscore = data_classify)
    }
    
    
    scorepar$n <- n+n_bg # recet after scoring
    
    return(scoresagainstclusters[,k])
    
    # }, mc.cores = k_clust)
  })
  
  for (kk in 1:k_clust){
    scoresagainstclusters[,kk] <- parRes[[kk]]
  }
  
  # assign cluster
  newallrelativeprobabsnotau <- graphClust:::allrelativeprobs(scoresagainstclusters)
  newallrelativeprobabs <- graphClust:::relativeprobswithtau(newallrelativeprobabsnotau,tauvec)
  
  newclustermembership<-graphClust:::reassignsamples(newallrelativeprobabs)
  
  return(list("clustermembership"=newclustermembership,"allrelativeprobabs"=newallrelativeprobabs, "shared_vars"=shared_vars, "classified_data"=myData))
}



reclass <- get_classification(cluster_res, mut_cov_data)
survdata$group <- reclass$clustermembership
survdata$Cluster <- survdata$group
survdata$Cluster <- as.factor(survdata$Cluster)


p <- c()

for(i in unique(survdata$ICC)){
  subdata <- survdata[survdata$ICC== i,]
  group_count <- subdata %>% 
    group_by(Cluster) %>% 
    dplyr::summarise(n = n())
  
  filtered_groups <- group_count %>% 
    filter(n >= 10) %>% 
    pull(Cluster)
  
  # Create a new DataFrame containing only the rows where group appears 5 or more times
  subdata <- subdata %>% 
    filter(Cluster %in% filtered_groups)
  
  # merge clinical information and cluster membership
  clinical <- list()
  clinical$event <- c(subdata$AMLt_STATUS,subdata$AMLt_STATUS)
  clinical$time <- c(subdata$AMLt_YEARS,subdata$AMLt_YEARS)
  clinical$group <- c(subdata$ICC,subdata$Cluster)
  clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)
  
  # THESE NOW MATCH
  length(clinical$group)
  length(clinical$time)
  
  table(clinical$group)
  
  # Kaplan-Meier curve for groups
  
  nb.cols <- length(unique(cluster_res$clustermembership))
  
  # 
  library(viridis)
  #colours_clusters <- c("blue", "#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")
  colours_clusters <- c("blue", viridis(6))
  #plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
  #legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')
  mypal=colours_clusters[c(1, c(seq(1:9)+1)[seq(1,6) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])]
  dimensione=length(mypal)-1
  
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
  library("survminer")
  file<-paste0("../comparison/Supp_",gsub("MDS/AML ", "MDSAML_",i),"_AMLonly.png")
  png(filename = file, width = 2300, height = 2000, res = 280)
  p[i] <- ggsurvplot(
    os,                     # survfit object with calculated statistics.
    data = clinical,             # data used to fit survival curves.
    #palette = colours_clusters,#colours_clusters, # personalized colours
    palette = alpha(mypal,c(1,rep(0.55, dimensione))), 
    # legend.labs = sort(as.character(unique(survdata$Cluster))), 
    risk.table = F,       # show risk table.
    pval = T,             # show p-value of log-rank test.
    conf.int = F,         # show confidence intervals for 
    # point estimates of survival curves.
    xlim = c(0,10),         # was c(0,14.4) 
    # survival estimates.
    xlab = "Time in years",   # customize X axis label.
    ylab = "Transformation probability",   # customize X axis label.
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.fontsize = 3,
    pval.size =4,
    title= paste0("AMLt - MDS only - ", gsub(pattern = "_", replacement = " ", i), " (n=",length(clinical$event)/2,")"),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
    
  )
  print(p[i])
  
  dev.off()
}

pdf(file = "../analysis/figures/Supplementary_MDSAML_AMLt_MDSonly.pdf", width = 16, height = 10)
ggarrange(p$`MDS/AML MR-Cyto`,p$`MDS/AML MR-Gen`, p$`MDS/AML TP53`, p$`MDS/AML NOS`, ncol = 2, nrow = 2 , labels = letters[1:4])
dev.off() 


