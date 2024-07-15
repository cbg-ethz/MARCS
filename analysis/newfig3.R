rm(list=ls())

library(graphClust)
library(scales)
library(ggplot2)
library(survival)
library(dplyr)
library(survival)
library(RColorBrewer)


cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")

table(cluster_res$clustermembership)


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

# calculate time to last follow-up
# survdata$time <- as.numeric(difftime(survdata$last_fu, survdata$firstdiagnosis_date, units="days"))
#survdata$time <- survdata$OS

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- survdata$OS_STATUS
clinical$time <- survdata$OS
clinical$group <- cluster_res$clustermembership
clinical$type <- survdata$WHO_2016

# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

# clinical$age <- #here needs to go the age of the patients#


# Change group from 1:22 to A:V
clinical$group <- as.factor(clinical$group)
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(clinical$group) <- lab#LETTERS[1:length(unique(clinical$group))]


table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))
colourysdots <- alpha(c("#202020","#771155","#114477","#ed2124","#771122","#DDDD77","#DDAA77","#117777","#CC99BB"), 0.7)
# 
# 
colours_clusters <- c("#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #c("#000000","#772266ff","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')

# mycolor2 <- c("#202020","#771122","#AA4455","#CC99BB","#774411","#AA7744",
#               "#DDAA77","#DD7788","#DDCC77","#77CCCC","#114477","#4477AA",
#               "#77AADD","#AA4488","#CC99BB","#DDCC77","#777711","#AAAA44","#117744","#44AA77",
#               "#88CCAA","#1122AA","#44AA00","#77CC77", "#AAAA44","#DDDD77","blue", "red", "yellow") ## COLOURS NEEDS TO BE ADJUSTED 


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

asd <- c("#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #c("#00000080","#771122","#11777780","#77112280","#11447780","#CC99BB80","#88CCAA80","#11774480","#77AADD80")

pos <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = asd, #colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
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
  title="Overall Survival (n=7480)",
  legend.labs = lab, # in legend of risk table,
  font.legend = list(size = 12),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)


library(graphClust)
library(scales)
library(ggplot2)


library(survival)
library(RColorBrewer)


#### PLOT FOR CANCER TYPES 
clinical <- list()
clinical$event <- survdata$OS_STATUS
clinical$time <- survdata$OS
clinical$IPSSR_ELN <- survdata$IPSSR_ELN
clinical$WHO_2016 <- survdata$WHO_2016
clinical$WHO_2022 <- survdata$WHO_2022
clinical$ICC <- survdata$ICC
clinical$group <-  cluster_res$clustermembership
clinical$group <- as.factor(clinical$group)
ab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
levels(clinical$group) <- lab
# levels(clinical$group) <- LETTERS[1:length(unique(clinical$group))]

# mycolor2 <- c("#202020","#771122","#AA4455","#CC99BB","#774411","#AA7744",
#               "#DDAA77","#DD7788","#DDCC77","#77CCCC","#114477","#4477AA",
#               "#77AADD","#AA4488","#CC99BB","#DDCC77","#777711","#AAAA44","#117744","#44AA77",
#               "#88CCAA","#1122AA","#44AA00","#77CC77", "#AAAA44","#DDDD77","blue", "red", "yellow")

# survival plot for each cancer type
clinical_aml <- as.data.frame(clinical)[clinical$IPSSR_ELN %in% c("ELN2017_adverse","ELN2017_intermediate","ELN2017_favorable", "s-AML"),]

group_count <- clinical_aml %>% 
  group_by(group) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 20 times
filtered_groups <- group_count %>% 
  filter(n >= 100) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_aml <- clinical_aml %>% 
  filter(group %in% filtered_groups)


os_aml <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = filtered_clinical_aml)
km_aml1 <- ggsurvplot(
  os_aml,                     # survfit object with calculated statistics.
  data = filtered_clinical_aml,             # data used to fit survival curves.
  palette = colours_clusters[levels(clinical$group) %in% unique(filtered_clinical_aml$group)], # personalized colours
  #legend.labs = c("UHR","HR1","NPM1","HR2","INT1","HR3"),
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
  filter(n >= 100) %>% 
  pull(group)

# Create a new DataFrame containing only the rows where group appears 5 or more times
filtered_clinical_aml <- clinical_aml %>% 
  filter(group %in% filtered_groups)

os_aml <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = filtered_clinical_aml)

km_aml2 <- ggsurvplot(
  os_aml,                     # survfit object with calculated statistics.
  data = filtered_clinical_aml,             # data used to fit survival curves.
  palette = colours_clusters[levels(clinical$group) %in% unique(filtered_clinical_aml$group)], # personalized colours
  #legend.labs = sort(as.character(unique(filtered_clinical_aml$group))), 
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

colours_clusters <- c("#2f0000","#BB3E03","#CA6702","#EE9B00","#0A9396","#005F73") #c("#000000","#771122aa","#CC99BBcc","#88CCAA","#117744","#77AADD")


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
  legend.labs = c("UHR","HR2","HR3", "INT2","LR1","LR2"), # in legend of risk table,
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table 
  
)



pdf(file = "../../forposter_fig2.pdf", width = 16.7, height = 6)
ggarrange( pos$plot, ppfs$plot,
           ncol = 2, nrow = 1 , labels = letters[3:6])
dev.off() 


pdf(file = "./figures/newfig3.pdf",  width = 15, height = 9)
ggarrange( pos$plot, ppfs$plot, pos$table, ppfs$table,
           ncol = 2, nrow = 2 , labels = c(letters[1:2], ""), heights = c(1,0.4) )
dev.off() 




# pdf(file = "./figures/SuppS5.pdf",  width = 16.7, height = 8)
# ggarrange( km_aml1$plot,km_aml2$plot, km_aml1$table,  km_aml2$table,
#            ncol = 2, nrow = 2 , labels = c(letters[1:2], "", ""), heights = c(1,0.4) )
# dev.off()

