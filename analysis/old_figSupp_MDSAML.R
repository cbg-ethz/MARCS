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

library(survival)
library(RColorBrewer)

# Load classification by mutation profile (cluster assignment) # CLUSTERS 
#cluster_res <- readRDS("../results/euler_memberships_final_9.rds")
# cluster_res <- readRDS("./results/euler_memberships_final_26.rds")
cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")

# remove this line in the future (after package update on euler)
# cluster_res$clustermembership <- cluster_res$clustermembership[[1]]
# cluster_res <- readRDS("../results/custer_res.rds")

table(cluster_res$clustermembership)

# mut_cov_data <- read.csv("../data/undivided_binary_matrix.csv")

# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")

# survdata$last_fu <- as.Date(survdata$last_fu, format = "%d.%m.%Y")
# survdata$firstdiagnosis_date <- as.Date(survdata$firstdiagnosis_date, format = "%d.%m.%Y")
# survdata$relapse_date <- as.Date(survdata$relapse_date, format = "%d.%m.%Y")
# survdata$type <- as.factor(survdata$)
survdata$WHO_2016 <- as.factor(survdata$WHO_2016)
survdata$WHO_2022 <- as.factor(survdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(survdata$IPSSR_ELN)
survdata$ICC <- as.factor(survdata$ICC)
survdata$Cluster <- cluster_res$clustermembership
survdata$Cluster <- as.factor(survdata$Cluster)
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
levels(survdata$Cluster) <- lab #LETTERS[1:length(unique(survdata$Cluster))]
survdata <- survdata[grep("MDS/AML",survdata$ICC),]



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
  colours_clusters <- c("blue", "#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")
  
  #plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
  #legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')
  mypal=colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])]
  dimensione=length(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])])-1
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
  library("survminer")
  file<-paste0("./figures/Supp_",gsub("MDS/AML ", "MDSAML_",i),".png")
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
    title= paste0("Overall Survival - ", gsub(pattern = "_", replacement = " ", i), " (n=",length(clinical$event)/2,")"),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
    
  )
  print(p[i])
  
  dev.off()
}

pdf(file = "./figures/Supplementary_MDS_AML_OS.pdf", width = 11, height = 10)
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
cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")

# remove this line in the future (after package update on euler)
# cluster_res$clustermembership <- cluster_res$clustermembership[[1]]
# cluster_res <- readRDS("../results/custer_res.rds")

table(cluster_res$clustermembership)

# mut_cov_data <- read.csv("../data/undivided_binary_matrix.csv")

# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")

# survdata$last_fu <- as.Date(survdata$last_fu, format = "%d.%m.%Y")
# survdata$firstdiagnosis_date <- as.Date(survdata$firstdiagnosis_date, format = "%d.%m.%Y")
# survdata$relapse_date <- as.Date(survdata$relapse_date, format = "%d.%m.%Y")
# survdata$type <- as.factor(survdata$)
survdata$WHO_2016 <- as.factor(survdata$WHO_2016)
survdata$WHO_2022 <- as.factor(survdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(survdata$IPSSR_ELN)
survdata$ELN2022_IPSSM <- as.factor(survdata$ELN2022_IPSSM)
survdata$Cluster <- cluster_res$clustermembership
survdata$Cluster <- as.factor(survdata$Cluster)
levels(survdata$Cluster) <- LETTERS[1:length(unique(survdata$Cluster))]
survdata <- survdata[!(survdata$ELN2022_IPSSM) %in% "IPSSM_NA",]
pfsdata <- read.delim("../data/IPSSM_df_clinical.tsv")
pfsdata <- pfsdata[, c("ID","AMLt_YEARS","AMLt_STATUS","LFS_YEARS","LFS_STATUS")]
survdata <- merge(survdata,pfsdata, by = "ID")
survdata$WHO_2016 <- as.factor(survdata$WHO_2016)
survdata$WHO_2022 <- as.factor(survdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(survdata$IPSSR_ELN)
survdata$ICC <- as.factor(survdata$ICC)
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
levels(survdata$Cluster) <- lab #LETTERS[1:length(unique(survdata$Cluster))]
survdata <- survdata[grep("MDS/AML",survdata$ICC),]


# survival plot for each cancer type
clinical_aml <- as.data.frame(clinical)[clinical$IPSSR_ELN %in% c("ELN2017_adverse","ELN2017_intermediate","ELN2017_favorable", "s-AML"),]

group_count <- clinical_aml %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 20 times
filtered_groups <- group_count %>% 
  filter(n >= 20) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_aml <- clinical_aml %>% 
  filter(group %in% filtered_groups)


os_aml <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = filtered_clinical_aml)
km_aml1 <- ggsurvplot(
  os_aml,                     # survfit object with calculated statistics.
  data = filtered_clinical_aml,             # data used to fit survival curves.
  palette = colours_clusters[which(levels(clinical$group) %in% sort(as.character(unique(filtered_clinical_aml$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_aml$group))), 
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
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("OS for AML (n=",dim(clinical_aml)[1],")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_aml1

# survival plot for each cancer type
clinical_aml <- as.data.frame(clinical)[!clinical$IPSSR_ELN %in% c("ELN2017_adverse","ELN2017_intermediate","ELN2017_favorable", "s-AML"),]


group_count <- clinical_aml %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 20 times
filtered_groups <- group_count %>% 
  filter(n >= 40) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_aml <- clinical_aml %>% 
  filter(group %in% filtered_groups)

os_aml <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = filtered_clinical_aml)

km_aml2 <- ggsurvplot(
  os_aml,                     # survfit object with calculated statistics.
  data = filtered_clinical_aml,             # data used to fit survival curves.
  palette = colours_clusters[which(levels(clinical$group) %in% sort(as.character(unique(filtered_clinical_aml$group))))], # personalized colours
  legend.labs = sort(as.character(unique(filtered_clinical_aml$group))), 
  risk.table = TRUE,       # show risk table.
  pval = T,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  pval.coord = c(0, 0.15),
  font.title= 16,
  title=paste0("OS for MDS (n=",dim(clinical_aml)[1],")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
); km_aml2



#PFS 

# Load classification by mutation profile (cluster assignment) # CLUSTERS 
cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")
mutation_covariate_data <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
id_clus <- data.frame(ID=mutation_covariate_data$ID, Cluster=cluster_res$clustermembership)

table(cluster_res$clustermembership)

# mut_cov_data <- read.csv("../data/undivided_binary_matrix.csv")

# import survival data
survdata <- read.delim("../data/IPSSM_df_clinical.tsv")
id_clus <- id_clus[id_clus$ID %in% survdata$ID,]
survdata <- merge(survdata, id_clus, by = "ID")
survdata$Cluster <- as.factor(survdata$Cluster)
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
levels(survdata$Cluster) <- lab#LETTERS[1:length(unique(survdata$Cluster))]

## remove groups with < 20 patients 
survdata <- survdata %>%
  group_by(Cluster) %>%
  filter(n() >= 40)

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- survdata$AMLt_STATUS
clinical$time <- survdata$AMLt_YEARS
clinical$group <- survdata$Cluster
#clinical$type <- survdata$WHO_2016

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

colours_clusters <- c("#000000","#771122aa","#CC99BBcc","#88CCAA","#117744","#77AADD")

pfs <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

ppfs <- ggsurvplot(
  pfs,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,15),         # present narrower X axis, but not affect
  # survival estimates.
  ylab = "Transformation probabiliy",
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  title="Time to AML transformation (n=2370)",
  font.title= 16,
  font.legend = list(size = 12),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table 
  
)





pdf(file = "newSupp_figAMLMDS_OS.pdf",  width = 16.7, height = 8)
ggarrange( km_aml1$plot,km_aml2$plot, km_aml1$table,  km_aml2$table,
           ncol = 2, nrow = 2 , labels = c(letters[1:2], "", ""), heights = c(1,0.4) )
dev.off()


