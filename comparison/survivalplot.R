rm(list=ls())

library(graphClust)
library(scales)
library(ggplot2)
library(survival)
library(dplyr)
library(survival)
library(viridis)
library(RColorBrewer)

# cluster_res <- readRDS("./results/euler_memberships_final_26.rds")
cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")
clusterAML <- readRDS("./AMLonly_6clus.rds")
clusterMDS <- readRDS("./MDSonly_6clus.rds")

table(cluster_res$clustermembership)
table(clusterAML$clustermembership)
table(clusterMDS$clustermembership)

AMLmat <- as.data.frame(clusterAML$clustermembership)
MDSmat <- as.data.frame(clusterMDS$clustermembership)
AMLmat$group <- paste0(AMLmat[,1],"_AML")
MDSmat$group <- paste0(MDSmat[,1],"_MDS")
colnames(AMLmat) <- c("cluster", "group")
colnames(MDSmat) <- c("cluster", "group")
combimembers <- rbind(AMLmat,MDSmat)
combimembers$ID <- rownames(combimembers)

TOTmat <- data.frame(ID=rownames(cluster_res$data), cluster=cluster_res$clustermembership, group=paste0(cluster_res$clustermembership))
TOTmat$group <- as.factor(TOTmat$group)
levels(TOTmat$group) <-LETTERS[1:9]
combimembers <- merge(TOTmat,combimembers, by = "ID")

# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k.csv")
totdata <- merge(survdata,combimembers, by="ID")

totdata$WHO_2016 <- as.factor(totdata$WHO_2016)
totdata$WHO_2022 <- as.factor(totdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(totdata$IPSSR_ELN)
totdata$ELN2022_IPSSM <- as.factor(survdata$ELN2022_IPSSM)
totdata$group.y <- as.factor(totdata$group.y)
totdata$group.y <- factor(totdata$group.y, levels = c("1_AML", "2_AML", "3_AML", "4_AML", "5_AML", "1_MDS", "2_MDS", "3_MDS", "4_MDS", "5_MDS", "6_MDS"))

# calculate time to last follow-up
# survdata$time <- as.numeric(difftime(survdata$last_fu, survdata$firstdiagnosis_date, units="days"))
#survdata$time <- survdata$OS

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- totdata$OS_STATUS
clinical$time <- totdata$OS
clinical$group <- totdata$group.y


# THESE NOW MATCH
length(clinical$group)
length(clinical$time)

# clinical$age <- #here needs to go the age of the patients#


# Change group from 1:22 to A:V
clinical$group <- as.factor(clinical$group)

#levels(clinical$group) <- LETTERS[1:length(unique(clinical$group))]

#table(clinical$group)

# Kaplan-Meier curve for groups

nb.cols <- length(unique(cluster_res$clustermembership))
colourysdots <- alpha(c("#202020","#771155","#114477","#ed2124","#771122","#DDDD77","#DDAA77","#117777","#CC99BB"), 0.7)
# 
# 
colours_clusters <- c("#000000","#772266ff","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")

plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')

mycolor2 <- c("#202020","#771122","#AA4455","#DD7788","#CC99BB",
             "#114477","#4477AA","#77AADD", "#77CCCC","#88CCAA","#117744")

# 
#               ,"#AA4488","#CC99BB","#DDCC77","#777711","#AAAA44",,"#44AA77",
#               ,"#1122AA","#44AA00","#77CC77", "#AAAA44","#DDDD77","blue", "red", "yellow") ## COLOURS NEEDS TO BE ADJUSTED 


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
library(viridis)
pos1 <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = viridis(11), #mycolor2,#colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,14.4),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title="Overall Survival (n=7480)",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)



# merge clinical information and cluster membership
clinical <- list()
clinical$event <- totdata$OS_STATUS
clinical$time <- totdata$OS
clinical$group <- combimembers$group.x
# THESE NOW MATCH
length(clinical$group)
length(clinical$time)
# Change group from 1:22 to A:V
clinical$group <- as.factor(clinical$group)

os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

pos2 <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = colours_clusters,#colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,14.4),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title="Overall Survival (n=7480)",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)




#
totmds <- totdata[grep("MDS",totdata$group.y),]
group_count <- totmds %>% 
  group_by(group.x) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 20 times
filtered_groups <- group_count %>% 
  filter(n >= 20) %>% 
  pull(group.x)

# Create a new DataFrame containing only the rows where group appears 5 or more times
totmds <- totmds %>% 
  filter(group.x %in% filtered_groups)



clinical <- list()
clinical$event <- totmds$OS_STATUS
clinical$time <- totmds$OS
clinical$group <- totmds$group.y


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

mdsos1 <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = viridis(6),#mycolor2[c(6:11)],#colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,14.4),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("Overall Survival - MDS only (n=", length(clinical$group), ")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)



# merge clinical information and cluster membership
clinical <- list()
clinical$event <- totmds$OS_STATUS
clinical$time <- totmds$OS
clinical$group <- totmds$group.x
# THESE NOW MATCH
length(clinical$group)
length(clinical$time)
# Change group from 1:22 to A:V
clinical$group <- as.factor(clinical$group)

os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

mdsos2 <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = colours_clusters[c(1,2,4,5,6,7,8,9)],#colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,14.4),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("Overall Survival - Combined (n=", length(clinical$group), ")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)








#
totaml <- totdata[grep("AML",totdata$group.y),]
group_count <- totaml %>% 
  group_by(group.x) %>% 
  dplyr::summarise(n = n())

# Filter out the groups that appear less than 20 times
filtered_groups <- group_count %>% 
  filter(n >= 20) %>% 
  pull(group.x)

# Create a new DataFrame containing only the rows where group appears 5 or more times
totaml <- totaml %>% 
  filter(group.x %in% filtered_groups)



clinical <- list()
clinical$event <- totaml$OS_STATUS
clinical$time <- totaml$OS
clinical$group <- totaml$group.y


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

amlos1 <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = plasma(5),#mycolor2[c(1:5)],#colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,14.4),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("Overall Survival - AML only (n=", length(clinical$group), ")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)



# merge clinical information and cluster membership
clinical <- list()
clinical$event <- totaml$OS_STATUS
clinical$time <- totaml$OS
clinical$group <- totaml$group.x
# THESE NOW MATCH
length(clinical$group)
length(clinical$time)
# Change group from 1:22 to A:V
clinical$group <- as.factor(clinical$group)

os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

amlos2 <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = colours_clusters[c(1,2,3,4,5,6,7)],#colours_clusters, # personalized colours
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,14.4),         # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("Overall Survival - Combined (n=", length(clinical$group), ")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)




#
# import survival data
pfsdata <- read.delim("../data/IPSSM_df_clinical.tsv")
pfsdata <- merge(totmds, pfsdata, by = "ID")

## remove groups with < 20 patients 
pfsdata <- pfsdata %>%
  group_by(group.x) %>%
  filter(n() >= 20)

# merge clinical information and cluster membership
clinical <- list()
clinical$event <- pfsdata$AMLt_STATUS
clinical$time <- pfsdata$AMLt_YEARS
clinical$group <- pfsdata$group.x
#clinical$type <- survdata$WHO_2016


pfs <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

pfs2 <- ggsurvplot(
  pfs,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = colours_clusters[c(1,4,6,7,8,9)], # personalized colours
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
  title=paste0("Time to AML transformation - Combined (n=", length(clinical$group), ")"),
  font.title= 16,
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table 
  
)


##
clinical <- list()
clinical$event <- pfsdata$AMLt_STATUS
clinical$time <- pfsdata$AMLt_YEARS
clinical$group <- pfsdata$group.y
clinical$groupSep <- pfsdata$group.y
clinical$groupTot <- pfsdata$group.x

#clinical$type <- survdata$WHO_2016


tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  groupSep , data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + groupSep , data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(unique(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("TRANSFORMATION Combined vs Separated Total: ", paste(round(LR, 1), "&", pvalue, sep = " ")))



pfs <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

pfs1 <- ggsurvplot(
  pfs,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = viridis(6),#mycolor2[c(6:11)], # personalized colours
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
  title=paste0("Time to AML transformation - MDS only (n=", length(clinical$group), ")"),
  font.title= 16,
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table 
  
)


pdf(file = "../figures/Supp_compare_OStot.pdf",  width = 16.7, height = 8)
ggarrange( pos2$plot, pos1$plot, pos2$table, pos1$table,
           ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
dev.off() 

pdf(file = "../figures/Supp_compare_OSAML.pdf",  width = 16.7, height = 8)
ggarrange( amlos2$plot, amlos1$plot, amlos2$table, amlos1$table,
           ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
dev.off() 

pdf(file = "../figures/Supp_compare_Tot.pdf",  width = 16.7, height = 22)
ggarrange( amlos2$plot, amlos1$plot, amlos2$table, amlos1$table,
           mdsos2$plot, mdsos1$plot, mdsos2$table, mdsos1$table,
           pfs2$plot, pfs1$plot, pfs2$table, pfs1$table,
           ncol = 2, nrow = 6 , labels = c(letters[1:2], "", "",letters[3:4], "", "",letters[5:6],"",""), heights = c(1,0.4,1,0.4,1,0.4) )
dev.off() 

pdf(file = "../figures/Supp_compare_OSMDS.pdf",  width = 16.7, height = 8)
ggarrange( mdsos2$plot, mdsos1$plot, mdsos2$table, mdsos1$table,
           ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
dev.off() 


pdf(file = "../figures/Supp_compare_PFS.pdf",  width = 16.7, height = 8)
ggarrange( pfs2$plot, pfs1$plot, pfs2$table, pfs1$table,
           ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
dev.off() 





pos1
pos2

amlos1
amlos2  

mdsos1
mdsos2

pfs1
pfs2


# ggarrange( pos1$plot, ppfs$plot, pos$table, ppfs$table,
#            km_aml1$plot,km_aml2$plot, km_aml1$table,  km_aml2$table,
#            ncol = 2, nrow = 4 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )


clinical <- list()
clinical$age <- totdata$AGE
clinical$sex <- totdata$SEX
clinical$event <- totdata$OS_STATUS
clinical$time <- totdata$OS
clinical$groupTot <- totdata$group.x
clinical$groupSep <- totdata$group.y
clinical$groupTot <- as.factor(clinical$groupTot)
clinical$groupSep <- as.factor(clinical$groupSep)

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + groupSep + age + sex, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("Combined vs Separated Total: ", paste(round(LR, 1), "&", pvalue, sep = " ")))



clinical <- list()
clinical$age <- totaml$AGE
clinical$sex <- totaml$SEX
clinical$event <- totaml$OS_STATUS
clinical$time <- totaml$OS
clinical$groupTot <- totaml$group.x
clinical$groupSep <- totaml$group.y
clinical$groupTot <- as.factor(clinical$groupTot)
clinical$groupSep <- as.factor(clinical$groupSep)

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + groupSep + age + sex, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("Combined vs Separated AML: ", paste(round(LR, 1), "&", pvalue, sep = " ")))


library(survival)
fitSep <- coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit")
fitTot <- coxph(Surv(time, as.numeric(event)) ~  groupTot+ age + sex, data = clinical, na.action = "na.omit")

concordance(fitSep)
concordance(fitTot)

# tdata <- as.data.frame(clinical) 
# tdata$time <- 10 # 5 years 
# pSep <- predict(fitSep, newdata =  tdata)
# pTot <- predict(fitTot)
# cfit <- concordance(Surv(time, as.numeric(event)) ~ pSep +  pTot , clinical)
# cfit
# round(coef(cfit), 3)
# round(cov2cor(vcov(cfit)), 3)  # high correlation




clinical <- list()
clinical$age <- totmds$AGE
clinical$sex <- totmds$SEX
clinical$event <- totmds$OS_STATUS
clinical$time <- totmds$OS
clinical$groupTot <- totmds$group.x
clinical$groupSep <- totmds$group.y
clinical$groupTot <- as.factor(clinical$groupTot)
clinical$groupSep <- as.factor(clinical$groupSep)

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + groupSep + age + sex, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("Combined vs Separated MDS: ", paste(round(LR, 1), "&", pvalue, sep = " ")))

library(survival)
fitSep <- coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit")
fitTot <- coxph(Surv(time, as.numeric(event)) ~  groupTot+ age + sex, data = clinical, na.action = "na.omit")

print("Separated MDS: ")
concordance(fitSep)
print("Combined MDS: ")
concordance(fitTot)
