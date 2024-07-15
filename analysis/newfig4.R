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
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(survdata$Cluster) <- lab #LETTERS[1:length(unique(survdata$Cluster))]
survdata <- survdata[!(survdata$ELN2022_IPSSM) %in% "IPSSM_NA",]

p <- c()

for(i in unique(survdata$ELN2022_IPSSM)){
  subdata <- survdata[survdata$ELN2022_IPSSM== i,]
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
  clinical$group <- c(subdata$ELN2022_IPSSM,subdata$Cluster)
  clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)
  
  # THESE NOW MATCH
  length(clinical$group)
  length(clinical$time)
  
  table(clinical$group)
  
  # Kaplan-Meier curve for groups
  
  nb.cols <- length(unique(cluster_res$clustermembership))
  
  # 
  # 
  colours_clusters <- c("red","#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") # "#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")
  
  #plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
  #legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')
  
  
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
  library("survminer")
  file<-paste0("../results/9_KM_comparison_",i,".png")
  png(filename = file, width = 2300, height = 2000, res = 280)
  p[i] <- ggsurvplot(
    os,                     # survfit object with calculated statistics.
    data = clinical,             # data used to fit survival curves.
    #palette = colours_clusters,#colours_clusters, # personalized colours
    palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
    # legend.labs = sort(as.character(unique(survdata$Cluster))), 
    risk.table = F,       # show risk table.
    legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
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
    pval.size =4.5,
    font.title= 14,
    pval.coord = c(0, 0.07),
    font.tickslab = c(13),
    font.legend = list(size = 10),
    title= paste0("OS - ", gsub(pattern = "_", replacement = " ", i)," (n=",length(clinical$event)/2,")"),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
    
  )
  print(p[i])
  
  dev.off()
}


#do ad hoc for ELN adverse
i <- "ELN2022_adverse"
subdata <- survdata[survdata$ELN2022_IPSSM== i,]
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
clinical$group <- c(subdata$ELN2022_IPSSM,subdata$Cluster)
clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))

# 
# 
colours_clusters <- c("red", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

#plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
#legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
file<-paste0("../results/9_KM_comparison_",i,".png")
png(filename = file, width = 2300, height = 2000, res = 280)
p[i] <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  #palette = colours_clusters,#colours_clusters, # personalized colours
  palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
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
  pval.coord = c(0, 0.07),
  pval.size =4.5,
  font.title= 14,
  font.legend = list(size = 10),
  legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
  font.tickslab = c(13),
  title= paste0("OS - ", gsub(pattern = "_", replacement = " ", i), " (n=",length(clinical$event)/2,")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)
print(p[i])

dev.off()





pdf(file = "../results/compare_strat_OS.pdf", width = 16, height = 16)
ggarrange(p$ELN2022_favorable, p$ELN2022_intermediate, p$ELN2022_adverse,
          p$`IPSSM_Very-Low`, p$IPSSM_Low, p$`IPSSM_Moderate-Low`,
          p$`IPSSM_Moderate-High`, p$IPSSM_High, p$`IPSSM_Very-High`, ncol = 3, nrow = 3 , labels = letters[1:9])
dev.off() 



### MERGE COUPLE OF RISK 
subdata <- survdata[survdata$ELN2022_IPSSM %in% c("IPSSM_Very-Low","IPSSM_Low"),]
group_count <- subdata %>% 
  group_by(Cluster) %>% 
  dplyr::summarise(n = n())

filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(Cluster)

# Create a new DataFrame containing only the rows where group appears 5 or more times
subdata <- subdata %>% 
  filter(Cluster %in% filtered_groups)

subdata$ELN2022_IPSSM <- as.factor("IPSSM_L_VL")

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- c(subdata$OS_STATUS,subdata$OS_STATUS)
clinical$time <- c(subdata$OS,subdata$OS)
clinical$group <- c(subdata$ELN2022_IPSSM, subdata$Cluster)
clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))

# 
# 
colours_clusters <- c("red","#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

#plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
#legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')

a <- c()

os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
file<-paste0("../results/9_KM_comparison_VL_L",".png")
png(filename = file, width = 2300, height = 2000, res = 280)
a[1] <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = "IPSSM_VL_L", invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
  #palette = alpha(c("red", "#CC99BBcc","#88CCAA","#117744","#77AADD"), c(rep(0.55, 4),1)), 
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
  legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
  pval.size =4.5,
  font.title= 14,
  pval.coord = c(0, 0.07),
  font.tickslab = c(13),
  font.legend = list(size = 10),
  title= paste0("OS - IPSSM Low, Very-Low", " (n=",length(clinical$event)/2,")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)
print(a[1])
dev.off()


### MERGE COUPLE OF RISK 
subdata <- survdata[survdata$ELN2022_IPSSM %in% c("IPSSM_Moderate-Low","IPSSM_Moderate-High"),]
group_count <- subdata %>% 
  group_by(Cluster) %>% 
  dplyr::summarise(n = n())

filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(Cluster)

# Create a new DataFrame containing only the rows where group appears 5 or more times
subdata <- subdata %>% 
  filter(Cluster %in% filtered_groups)

subdata$ELN2022_IPSSM <- as.factor("IPSSM_MH_ML")

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- c(subdata$OS_STATUS,subdata$OS_STATUS)
clinical$time <- c(subdata$OS,subdata$OS)
clinical$group <- c(subdata$ELN2022_IPSSM, subdata$Cluster)
clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))

# 
# 
colours_clusters <- c("red", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")#"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

#plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
#legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
file<-paste0("../results/9_KM_comparison_ML_MH",".png")
png(filename = file, width = 2300, height = 2000, res = 280)
a[2] <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = "IPSSM_ML_MH", invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
  #palette = alpha(c("red", "#CC99BBcc","#88CCAA","#117744","#77AADD"), c(rep(0.55, 4),1)), 
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
  pval.size =4.5,
  legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
  font.title= 14,
  pval.coord = c(0, 0.07),
  font.tickslab = c(13),
  font.legend = list(size = 10),
  title= paste0("OS - IPSSM  Moderate-High, Moderate-Low", " (n=",length(clinical$event)/2,")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)
print(a[2])
dev.off()



### MERGE COUPLE OF RISK 
subdata <- survdata[survdata$ELN2022_IPSSM %in% c("IPSSM_Very-High","IPSSM_High"),]
group_count <- subdata %>% 
  group_by(Cluster) %>% 
  dplyr::summarise(n = n())

filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(Cluster)

# Create a new DataFrame containing only the rows where group appears 5 or more times
subdata <- subdata %>% 
  filter(Cluster %in% filtered_groups)

subdata$ELN2022_IPSSM <- as.factor("IPSSM_VH_H")

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- c(subdata$OS_STATUS,subdata$OS_STATUS)
clinical$time <- c(subdata$OS,subdata$OS)
clinical$group <- c(subdata$ELN2022_IPSSM, subdata$Cluster)
clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))

# 
# 
colours_clusters <- c("red", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

#plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
#legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
file<-paste0("../results/9_KM_comparison_H_VH",".png")
png(filename = file, width = 2300, height = 2000, res = 280)
a[3] <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = "IPSSM_H_VH", invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
  #palette = alpha(c("red", "#CC99BBcc","#88CCAA","#117744","#77AADD"), c(rep(0.55, 4),1)), 
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
  pval.size =4.5,
  font.title= 14,
  font.tickslab = c(13),
  legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
  pval.coord = c(0, 0.07),
  font.legend = list(size = 10),
  title= paste0("OS - IPSSM Very-High, High", " (n=",length(clinical$event)/2,")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)
print(a[3])
dev.off()

pdf(file = "./figures/newfig4.pdf", width = 16.7, height = 12)
ggarrange( p$ELN2022_adverse, 
           p$ELN2022_intermediate,p$ELN2022_favorable, a[[3]],a[[2]],a[[1]], 
           ncol = 3, nrow = 2 , labels = letters[1:6])
dev.off() 


pdf(file = "../../forposter_fig3.pdf", width = 16.7, height = 5)
ggarrange( p$ELN2022_adverse, 
           p$ELN2022_intermediate,p$ELN2022_favorable,
           ncol = 3, nrow = 1 , labels = letters[3:6])
dev.off() 

##### CREATE A COMBI FIGURE

library(graphClust)
library(scales)
library(ggplot2)
library(survival)
library(dplyr)
# Load required libraries
library(survival)
library(survminer)


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
survdata$ELN2022_IPSSM <- as.factor(survdata$ELN2022_IPSSM)
survdata$Cluster <- cluster_res$clustermembership
survdata$Cluster <- as.factor(survdata$Cluster)
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
levels(survdata$Cluster) <- lab #LETTERS[1:length(unique(survdata$Cluster))]
survdata <- survdata[!(survdata$ELN2022_IPSSM) %in% "IPSSM_NA",]
pfsdata <- read.delim("../data/IPSSM_df_clinical.tsv")
pfsdata <- pfsdata[, c("ID","AMLt_YEARS","AMLt_STATUS","LFS_YEARS","LFS_STATUS")]
survdata <- merge(survdata,pfsdata, by = "ID")


### MERGE COUPLE OF RISK 
subdata <- survdata[survdata$ELN2022_IPSSM %in% c("IPSSM_Very-Low","IPSSM_Low"),]
group_count <- subdata %>% 
  group_by(Cluster) %>% 
  dplyr::summarise(n = n())

filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(Cluster)

# Create a new DataFrame containing only the rows where group appears 5 or more times
subdata <- subdata %>% 
  filter(Cluster %in% filtered_groups)

subdata$ELN2022_IPSSM <- as.factor("IPSSM_VL_L")

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- c(subdata$AMLt_STATUS,subdata$AMLt_STATUS)
clinical$time <- c(subdata$AMLt_YEARS,subdata$AMLt_YEARS)
clinical$group <- c(subdata$ELN2022_IPSSM, subdata$Cluster)
clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))

# 
# 
colours_clusters <- c("red", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")#"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

b <- c()

os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
file<-paste0("../results/9_KM_tAML_comparison_VL_L",".png")
png(filename = file, width = 2300, height = 2000, res = 280)
b[1] <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = "IPSSM_VL_L", invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
  #palette = alpha(c("red", "#CC99BBcc","#88CCAA","#117744","#77AADD"), c(rep(0.55, 4),1)), 
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
  pval.size =4.5,   font.title= 14,
  font.legend = list(size = 10),
  pval.coord = c(0, 0.07),
  legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
  font.tickslab = c(11),
  title= paste0("AMLt - IPSSM Low, Very-Low", " (n=",length(clinical$event)/2,")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)
print(b[1])
dev.off()


### MERGE COUPLE OF RISK 
subdata <- survdata[survdata$ELN2022_IPSSM %in% c("IPSSM_Moderate-Low","IPSSM_Moderate-High"),]
group_count <- subdata %>% 
  group_by(Cluster) %>% 
  dplyr::summarise(n = n())

filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(Cluster)

# Create a new DataFrame containing only the rows where group appears 5 or more times
subdata <- subdata %>% 
  filter(Cluster %in% filtered_groups)

subdata$ELN2022_IPSSM <- as.factor("IPSSM_ML_MH")

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- c(subdata$AMLt_STATUS,subdata$AMLt_STATUS)
clinical$time <- c(subdata$AMLt_YEARS,subdata$AMLt_YEARS)
clinical$group <- c(subdata$ELN2022_IPSSM, subdata$Cluster)
clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))

# 
# 
colours_clusters <- c("red", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

#plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
#legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
file<-paste0("../results/9_KM_tAML_comparison_ML_MH",".png")
png(filename = file, width = 2300, height = 2000, res = 280)
b[2] <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = "IPSSM_ML_MH", invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
  #palette = alpha(c("red", "#CC99BBcc","#88CCAA","#117744","#77AADD"), c(rep(0.55, 4),1)), 
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
  pval.size =4.5,   font.title= 14,
  legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
  font.tickslab = c(11),
  pval.coord = c(0, 0.07),
  font.legend = list(size = 10),
  title= paste0("AMLt - IPSSM Moderate-High, Moderate-Low", " (n=",length(clinical$event)/2,")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)
print(b[2])
dev.off()



### MERGE COUPLE OF RISK 
subdata <- survdata[survdata$ELN2022_IPSSM %in% c("IPSSM_Very-High","IPSSM_High"),]
group_count <- subdata %>% 
  group_by(Cluster) %>% 
  dplyr::summarise(n = n())

filtered_groups <- group_count %>% 
  filter(n >= 10) %>% 
  pull(Cluster)

# Create a new DataFrame containing only the rows where group appears 5 or more times
subdata <- subdata %>% 
  filter(Cluster %in% filtered_groups)

subdata$ELN2022_IPSSM <- as.factor("IPSSM_H_VH")

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- c(subdata$OS_STATUS,subdata$OS_STATUS)
clinical$time <- c(subdata$OS,subdata$OS)
clinical$group <- c(subdata$ELN2022_IPSSM, subdata$Cluster)
clinical$type <- c(subdata$WHO_2016,subdata$WHO_2016)

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))

# 
# 
colours_clusters <- c("red", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

#plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
#legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
file<-paste0("../results/9_KM_tAML_comparison_H_VH",".png")
png(filename = file, width = 2300, height = 2000, res = 280)
b[3] <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = alpha(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = "IPSSM_H_VH", invert = T)]])],c(1,rep(0.55, length(unique(subdata$Cluster))))), 
  #palette = alpha(c("red", "#CC99BBcc","#88CCAA","#117744","#77AADD"), c(rep(0.55, 4),1)), 
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
  ylab = "Transformation probability",   # customize X axis label.
  pval.size =4.5,  font.title= 14,
  font.tickslab = c(11),
  pval.coord = c(0, 0.07),
  font.legend = list(size = 10),
  legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
  title= paste0("AMLt - IPSSM Very-High, High", " (n=",length(clinical$event)/2,")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)
print(b[3])
dev.off()


pdf(file = "./figures/combifig4.pdf", width = 16.7, height = 18)
ggarrange( p$ELN2022_adverse, 
           p$ELN2022_intermediate,p$ELN2022_favorable, a[[3]],a[[2]],a[[1]],
           b[[3]],b[[2]],b[[1]],
           ncol = 3, nrow = 3 , labels = letters[1:9])
dev.off() 

pdf(file = "../../forposter_fig_combi.pdf", width = 16.7, height = 10)
ggarrange( p$ELN2022_adverse, 
           p$ELN2022_intermediate,p$ELN2022_favorable,a[[3]], b[[3]],b[[1]],
           ncol = 3, nrow = 2 , labels = letters[1:9])
dev.off() 
