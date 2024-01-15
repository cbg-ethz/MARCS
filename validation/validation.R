library(clustNet)
library(parallel)
library(xlsx)

rm(list=ls())

cluster_results <- readRDS("../results/euler_memberships_8k_9clusters.rds")


validata <- read.table("./mda_binary-mutationCovariate-matrix.txt")[,-58]
mutation_covariate_data <- readRDS("../validation//aml_data.rds")
validata$UPI <- mutation_covariate_data$UPI
riskstrat <- readRDS("../validation//risk_scores.rds")

# ids <- mutation_covariate_data$UPI[mutation_covariate_data$Dx %in% c("AML","MDS")]
# validata <- validata[validata$UPI %in% ids,]

# totdata <- merge(totdata, riskstrat[,c(1,85)], by="UPI")
# totdata <- totdata[totdata$UPI %in% ids,]

colnames(validata)[colnames(validata)=="CG_8"] <- "X.8"
colnames(validata)[colnames(validata)=="CG_5"] <- "X.5"
colnames(validata)[colnames(validata)=="CG_7"] <- "X.7"
colnames(validata)[colnames(validata)=="CG_inv16"] <- "inv.16."
colnames(validata)[colnames(validata)=="CG_8_21"] <- "t.8.21."
colnames(validata)[colnames(validata)=="BM_Blast"] <- "BM_BLASTS"
colnames(validata)[colnames(validata)=="PB_HGB"] <- "HB"
colnames(validata)[colnames(validata)=="PB_PLT"] <- "PLT"
colnames(validata)[colnames(validata)=="PB_WBC"] <- "WBC"
colnames(validata)[colnames(validata)=="age"] <- "AGE"
colnames(validata)[colnames(validata)=="sex"] <- "SEX"

dim(validata)


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

rownames(validata) <- validata$UPI

mds_reclass <- get_classification(cluster_results, validata)

### GENE BARPLOT FOR RECLASSIFIED DATA 
mutation_select <-validata
sharedids <- colnames(mutation_select)[colnames(mutation_select) %in% colnames(cluster_results$data)]
subdata <- mutation_select[,sharedids]
subdata$SEX <- NULL 
subdata$AGE <- NULL 
subdata$BM_BLASTS <- NULL 
subdata$PLT <- NULL 
subdata$WBC <- NULL 
subdata$HB <- NULL 
colnames(subdata)[32:36]<- c("-5","-7","+8","inv(16)", "t(8;21)")
subdata <- cbind(subdata,mds_reclass$clustermembership)
mutation_sum <- aggregate(subdata, by = list(subdata$`mds_reclass$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`mds_reclass$clustermembership` <- NULL
# Reshape the data for plotting
library(ggplot2)
library(reshape2)
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
levels(melted_data_sum$Cluster) <- LETTERS[1:length(unique(cluster_results$clustermembership))]

melted_data_sum$Cluster[melted_data_sum$Cluster=="1"] <- "A"
melted_data_sum$Cluster[melted_data_sum$Cluster=="2"] <- "B"
melted_data_sum$Cluster[melted_data_sum$Cluster=="3"] <- "C"
melted_data_sum$Cluster[melted_data_sum$Cluster=="4"] <- "D"
melted_data_sum$Cluster[melted_data_sum$Cluster=="5"] <- "E"
melted_data_sum$Cluster[melted_data_sum$Cluster=="6"] <- "F"
melted_data_sum$Cluster[melted_data_sum$Cluster=="7"] <- "G"
melted_data_sum$Cluster[melted_data_sum$Cluster=="8"] <- "H"
melted_data_sum$Cluster[melted_data_sum$Cluster=="9"] <- "I"


colours_clusters <- c("#000000","#772266","#117777","#7c1a29","#114477","#cd9bbc","#88CCAA","#117744","#77AADD")
# Stacked barplot (sum of mutations)
barplot_sum0 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("A"="#000000","B"="#772266","C"="#117777","D"="#7c1a29","E"="#114477","F"="#cd9bbc","G"= "#88CCAA","H"="#117744","I"="#77AADD"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =12))+
  theme(axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 14))+
  ggtitle(label = paste("Genetics")) + theme(plot.title = element_text(color="black", size=15, face="bold"))+
  theme(legend.text=element_text(size=10), legend.title = element_text(size = 12))



#### HERE DO THE KM PLOT 
library(survival)
library(RColorBrewer)
library(survminer)
library(dplyr)
library(ggpubr)
library(reshape2)

clinical <- readRDS("../validation/os_data.rds")

# read data
cluster_results <- mds_reclass
mutation_covariate_data <- readRDS("../validation//aml_data.rds")

# merge features
clinical$group <- as.factor(cluster_results$clustermembership)
levels(clinical$group) <- LETTERS[1:9]

clinical$type <- mutation_covariate_data$Dx
clinical$gender <- mutation_covariate_data$Gender
clinical$age <- mutation_covariate_data$age



# Kaplan-Meier curve for groups
colourysdots <- c("#202020","#774411","#DDAA77","#ed2124","#114477","#CC99BB",
                  "#88CCAA","#117744","#77AADD")


# p_surv_temp <- ggsurvplot(survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical), clinical, palette = colourysdots,legend = "none", xlab="Time (years)", size=0.6)
# 
# p_surv <- p_surv_temp$plot +scale_y_continuous(
#   limits = c(0, 1),
#   breaks = c(0,0.2,0.4,0.6,0.8,1),
#   expand = c(0.01, 0)
# );p_surv
# 
# # save plot
# saveRDS(p_surv, "../figures/km_plot.rds")




# mycolor2 <- c("#202020","#771122","#AA4455","#CC99BB","#774411","#AA7744",
#                        "#DDAA77","#DD7788","#DDCC77","#77CCCC","#114477","#4477AA",
#                        "#77AADD","#AA4488","#CC99BB","#DDCC77","#777711","#AAAA44","#117744","#44AA77",
#                        "#88CCAA","#1122AA","#44AA00","#77CC77", "#AAAA44","#DDDD77") ## COLOURS NEEDS TO BE ADJUSTED

# investigate subset of patients diagnosed with AML and MDS 
risk_scores <- readRDS("../validation//risk_scores.rds")

mean(risk_scores$PB_HGB[risk_scores$Dx=="MDS"])
mean(risk_scores$PB_HGB[risk_scores$Dx=="AML"])

clinical$risk_scores <- risk_scores$risk_scores
clinical$dx <- risk_scores$Dx 

clinical_aml_mds <- as.data.frame(clinical)
clinical_aml_mds <- clinical_aml_mds[clinical$dx %in% c("AML","MDS"),] #clinical_aml_mds[!is.na(clinical$risk_scores),]

# check which classification is more predictive in survival
# check IPSSR_ELN vs our clustering, given the WHO_2016
summary(coxph(Surv(time, as.numeric(event)) ~ group + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]
summary(coxph(Surv(time, as.numeric(event)) ~ risk_scores + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]

# Clinical + tissue
tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~ risk_scores + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]
# Clinical + tissue + group
groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ group + risk_scores + type + age + gender, data = clinical_aml_mds, na.action = "na.omit"))$logtest[1]

LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$group)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
paste(round(LR, 1), "&", pvalue, sep = " ")

colours_clusters <- c("#000000","#772266ff","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

os_aml_mds <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = clinical_aml_mds)

pos <- ggsurvplot(
  os_aml_mds,                     # survfit object with calculated statistics.
  data = clinical_aml_mds,             # data used to fit survival curves.
  palette = colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = T,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  # xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  font.title= 16,
  font.x=14,
  font.y=14,
  font.xtickslab=12,
  font.ytickslab=12,
  xlim = c(0,10),
  pval.size =4,
  legend.title="", 
  title="Overall survival on validation cohort (n=1035)",
  risk.table.y.text = FALSE # show bars instead of names in text annotations,
  # in legend of risk table
)



riskcolor <- c("#202020","#24878E","#88CCAA",
               "#453581","#DDAA77","#34618D",
               "#774411", "#440154", "#77CC77")



os_aml_mds <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ risk_scores, data = clinical_aml_mds)
pros <- ggsurvplot(
  os_aml_mds,                     # survfit object with calculated statistics.
  data = clinical_aml_mds,             # data used to fit survival curves.
  palette = riskcolor, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,10),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  font.title= 16,
  font.x=14,
  font.y=14,
  font.xtickslab=12,
  font.ytickslab=12,
  risk.table.fontsize = 4,
  pval.size =4,
  title="Corresponding IPSS-M/ELN2022 scores",
  risk.table.y.text = F, # show bars instead of names in text annotations
  # in legend of risk table
  legend.title="", 
  legend.labs=c(unique(clinical_aml_mds$risk_scores))
)




# survival plot for each cancer type
clinical_aml <- as.data.frame(clinical)[clinical$type=="AML",]

group_count <- clinical_aml %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 5 times
filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_aml <- clinical_aml %>% 
  filter(group %in% filtered_groups)

os_aml <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = filtered_clinical_aml)
km_aml <- ggsurvplot(
  os_aml,                     # survfit object with calculated statistics.
  data = filtered_clinical_aml,             # data used to fit survival curves.
  palette = colours_clusters[which(levels(clinical_aml_mds$group) %in% sort(as.character(unique(filtered_clinical_aml$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_aml$group))), 
  risk.table = TRUE,       # show risk table.
  pval = T,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,10),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  font.title= 16,
  font.x=14,
  font.y=14,
  font.xtickslab=12,
  font.ytickslab=12,
  risk.table.fontsize = 4,
  pval.size =4,
  title=paste0("OS for AML (n=",dim(filtered_clinical_aml)[1],")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_aml


clinical_mds <- as.data.frame(clinical)[clinical$type=="MDS",]

group_count <- clinical_mds %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 5 times
filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_mds <- clinical_mds %>% 
  filter(group %in% filtered_groups)

os_mds <- survfit(Surv(time = as.numeric(time)/12, event = as.numeric(event)) ~ group, data = filtered_clinical_mds)
km_mds <- ggsurvplot(
  os_mds,                     # survfit object with calculated statistics.
  data = filtered_clinical_mds,             # data used to fit survival curves.
  palette = colours_clusters[which(levels(clinical_aml_mds$group) %in% sort(as.character(unique(filtered_clinical_mds$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_mds$group))), 
  risk.table = TRUE,       # show risk table.
  pval = T,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,10),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  font.title= 16,
  font.x=14,
  font.y=14,
  font.xtickslab=12,
  font.ytickslab=12,
  risk.table.fontsize = 4,
  pval.size =4,
  title=paste0("OS for MDS (n=",dim(filtered_clinical_mds)[1],")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_mds








pdf(file = "../validation/validation.pdf", width = 16.7, height = 20)
ggarrange(ggarrange(pos$plot, pros$plot,
          pos$table, pros$table,
          km_aml$plot,km_mds$plot,
          km_aml$table,km_mds$table,
          nrow=4, ncol=2, labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(3,1,3,1)),
          barplot_sum0,
          ncol = 1, nrow = 2, 
          labels = c("","e"),
          heights = c(4,1))
dev.off()
