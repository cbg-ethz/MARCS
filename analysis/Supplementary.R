library(graphClust)
library(scales)
library(ggplot2)
library(ggpubr)
library(survival)
rm(list=ls())
library(survival)
library(RColorBrewer)

# Load classification by mutation profile (cluster assignment) # CLUSTERS 
#cluster_res <- readRDS("../results/euler_memberships_final_9.rds")
# cluster_res <- readRDS("./results/euler_memberships_final_26.rds")
cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")
table(cluster_res$clustermembership)

groups <- data.frame(cluster_res$clustermembership)
groups$ID <- rownames(cluster_res$data)
colnames(groups)[1] <- "Cluster"
groups$Cluster <- as.factor(groups$Cluster)
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(groups$Cluster) <- lab # LETTERS[1:length(unique(groups$Cluster))]

# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
colnames(survdata)[34:46] <- c("+8", "complex", "-17", "-5", "-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
data <- merge(survdata,groups,by="ID")

data$Cluster <- factor(data$Cluster, levels= c("UHR","HR1", "HR2", "HR3", "INT1", "INT2", "NPM1","LR1" , "LR2") )

library(dplyr)
library(ggplot2)

# Calculate the total number of mutations per patient
data <- data %>%
  mutate(total_mutations = rowSums(select(., ASXL1:NRAS, +8:"t(v;11)")))

# Calculate the average number of mutations per "Cluster"
avg_mutations_per_cluster <- data %>%
  group_by(Cluster) %>%
  summarize(avg_mutations = mean(total_mutations))

# View the average mutations per cluster
print(avg_mutations_per_cluster)

# Create a boxplot showing the number of mutations per "Cluster"
tmb <- ggplot(data, aes(x = Cluster, y = total_mutations, fill=Cluster)) +
  geom_boxplot() +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  labs(title = "Number of mutations per Cluster", 
       x = "Cluster", 
       y = "Number of mutations")+
  theme_minimal()




library(dplyr)

# Assuming data is the name of your data frame

avg_blasts_per_cluster <- data %>%
  group_by(Cluster) %>%
  summarize(avg_blasts = median(BM_BLASTS, na.rm = TRUE))

# Print the average BM_BLASTS per cluster
print(avg_blasts_per_cluster)

blast <- ggplot(data, aes(x = Cluster, y = BM_BLASTS, fill = Cluster)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  labs(title = "BM Blasts per Cluster", 
       x = "Cluster", 
       y = "Blasts (%)") +
  theme_minimal()
# png(filename = "../results/9_blasts.png", width = 1200, height = 700, res = 150)
# blast
# dev.off()
png(filename = "./figures/SuppS4.png", width = 3000, height = 1300, res = 200)
ggarrange(blast,tmb, labels = c("a","b"))
dev.off()

###### SUPPLEMENTARY S5 

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

os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

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
ab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(clinical$group) <- lab
# levels(clinical$group) <- LETTERS[1:length(unique(clinical$group))]

mycolor2 <- c("#202020","#771122","#AA4455","#CC99BB","#774411","#AA7744",
              "#DDAA77","#DD7788","#DDCC77","#77CCCC","#114477","#4477AA",
              "#77AADD","#AA4488","#CC99BB","#DDCC77","#777711","#AAAA44","#117744","#44AA77",
              "#88CCAA","#1122AA","#44AA00","#77CC77", "#AAAA44","#DDDD77","blue", "red", "yellow")

# survival plot for each cancer type
aml_vector <- c("AML t(8;21)", "AML NOS", "AML MRC", "AML NPM1", "AML CEBPAbi",
                "AML RUNX1", "AML inv(3)", "AML inv(16)", "AML t(6;9)",
                "AML t(9;11)", "AML t(15;17)", "aCML")
clinical_aml <- as.data.frame(clinical)[clinical$WHO_2016 %in% aml_vector,] #c("ELN2017_adverse","ELN2017_intermediate","ELN2017_favorable", "s-AML"),]

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
#clinical_aml <- as.data.frame(clinical)[!clinical$IPSSR_ELN %in% c("ELN2017_adverse","ELN2017_intermediate","ELN2017_favorable", "s-AML"),]
clinical_aml <- as.data.frame(clinical)[!clinical$WHO_2016 %in% aml_vector,] 

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

pdf(file = "./figures/SuppS5.pdf",  width = 16.7, height = 8)
ggarrange( km_aml1$plot,km_aml2$plot, km_aml1$table,  km_aml2$table,
           ncol = 2, nrow = 2 , labels = c(letters[1:2], "", ""), heights = c(1,0.4) )
dev.off()

### SUPPLEMENTARY S6 

# get bar plots of cluster assignment
rm(list=ls())

library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# read data
#cluster_results <- readRDS("../results/euler_memberships_final_9.rds")
cluster_results <- readRDS("../results/euler_memberships_8k_9clusters.rds")
table(cluster_results$clustermembership)

mutation_covariate_data <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
table(cluster_results$clustermembership)

#colnames(mutation_covariate_data)[43:65] <- c("+11","+13","+21", "+22", "+8", "complex", "-12", "-17","-18", "-20", "-3", "-5/-5q", "-7/-7q", "-9", "inv(16)","inv(3)",
#                                              "-Y", "t(15;17)", "t(6;9)", "t(8;21)",  "t(9;11)",  "t(9;22)", "t(v;11)")

colnames(mutation_covariate_data)[34:46] <- c("+8", "complex", "-17", "-5", "-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")

#colnames(mutation_covariate_data)[c(71:74,77,78)] <- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")

#colnames(cluster_results$data)[42:64] <- c("+11","+13","+21", "+22", "+8", "complex", "-12", "-17","-18", "-20", "-3", "-5/-5q", "-7/-7q", "-9", "inv(16)","inv(3)",
# "-Y", "t(15;17)", "t(6;9)", "t(8;21)",  "t(9;11)",  "t(9;22)", "t(v;11)")
#colnames(cluster_results$data)[c(65:70)] <- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")

colnames(cluster_results$data)[33:45] <- c("+8", "complex", "-17", "-5", "-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")


for(i in 1:length(cluster_results$DAGs)){
  rownames(cluster_results$DAGs[[i]])[33:45] <- c("+8", "complex", "-17", "-5", "-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
  #rownames(cluster_results$DAGs[[i]])[c(65:70)]<- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")
  
  colnames(cluster_results$DAGs[[i]])[33:45] <- c("+8", "complex", "-17", "-5", "-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
  #colnames(cluster_results$DAGs[[i]])[c(65:70)]<- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")
}





#FILL OF MISSING IPSS-R AND IPSS-M
#mutation_covariate_data$IPSSR_ELN[mutation_covariate_data$IPSSR_ELN=="IPSSR_NA"] <- "IPSSR_Int"
mutation_covariate_data$ELN2022_IPSSM[mutation_covariate_data$ELN2022_IPSSM=="IPSSM_NA"] <- "IPSSM_Moderate-High"

### GENE BARPLOTS 

# BARPLOT WITH MUTATIONS 
#FILL OF MISSING IPSS-R AND IPSS-M
library(ggplot2)
library(reshape2)
#mutation_covariate_data$IPSSR_ELN[mutation_covariate_data$IPSSR_ELN=="IPSSR_NA"] <- "IPSSR_Int"
mutation_covariate_data$ELN2022_IPSSM[mutation_covariate_data$ELN2022_IPSSM=="IPSSM_NA"] <- "IPSSM_Moderate-High"
mutation_covariate_data <- cbind(mutation_covariate_data,cluster_results$clustermembership)

mutation_select <- mutation_covariate_data
subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL

#mutation_covariate_data$IPSSR_ELN[mutation_covariate_data$IPSSR_ELN=="IPSSR_NA"] <- "IPSSR_Int"
mutation_covariate_data$ELN2022_IPSSM[mutation_covariate_data$ELN2022_IPSSM=="IPSSM_NA"] <- "IPSSM_Moderate-High"
mutation_covariate_data <- cbind(mutation_covariate_data,cluster_results$clustermembership)

mutation_select <- mutation_covariate_data
subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL
# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
levels(melted_data_sum$Cluster) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

colours_clusters <- c("#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #c("#000000","#772266","#117777","#7c1a29","#114477","#cd9bbc","#88CCAA","#117744","#77AADD")
# Stacked barplot (sum of mutations)


###1
classes <- c("IPSSM_Very-Low","IPSSM_Low") #c("IPSSM_Moderate-High","IPSSM_Moderate-Low")  #c("IPSSM_Very-High", "IPSSM_High")
mutation_select <- mutation_covariate_data[mutation_covariate_data$ELN2022_IPSSM %in% classes,]

subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL


# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
melted_data_sum$Cluster <- as.factor(melted_data_sum$Cluster)
levels(melted_data_sum$Cluster) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")



# Stacked barplot (sum of mutations)
barplot_sum1 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 10))+
  ggtitle(label = paste("Genetics for", classes, collapse = " ")) + theme(plot.title = element_text(color="black", size=10, face="bold"))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size = 10))

#png(filename = paste0("../results/9_genebarplot_",classes[1],".png"), width = 900, height = 500)
print(barplot_sum1)
#dev.off()

###2
classes <- c("IPSSM_Moderate-High","IPSSM_Moderate-Low")  #c("IPSSM_Very-High", "IPSSM_High")
mutation_select <- mutation_covariate_data[mutation_covariate_data$ELN2022_IPSSM %in% classes,]

subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL
# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
melted_data_sum$Cluster <- as.factor(melted_data_sum$Cluster)
levels(melted_data_sum$Cluster) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

# Stacked barplot (sum of mutations)
barplot_sum2 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 10))+
  ggtitle(label = paste("Genetics for", classes, collapse = " ")) + theme(plot.title = element_text(color="black", size=10, face="bold"))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size = 10))

# png(filename = paste0("../results/9_genebarplot_",classes[1],".png"), width = 900, height = 500)
print(barplot_sum2)
# dev.off()


###3
classes <- c("IPSSM_Very-High", "IPSSM_High")
mutation_select <- mutation_covariate_data[mutation_covariate_data$ELN2022_IPSSM %in% classes,]

subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL
# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
melted_data_sum$Cluster <- as.factor(melted_data_sum$Cluster)
levels(melted_data_sum$Cluster) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

# Stacked barplot (sum of mutations)
barplot_sum3 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 10))+
  ggtitle(label = paste("Genetics for", classes, collapse = " ")) + theme(plot.title = element_text(color="black", size=10, face="bold"))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size = 10))

# png(filename = paste0("../results/9_genebarplot_",classes[1],".png"), width = 900, height = 500)
print(barplot_sum3)
# dev.off()

###4
classes <- c("ELN2022_adverse")
mutation_select <- mutation_covariate_data[mutation_covariate_data$ELN2022_IPSSM %in% classes,]

subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL
# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
melted_data_sum$Cluster <- as.factor(melted_data_sum$Cluster)
levels(melted_data_sum$Cluster) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

# Stacked barplot (sum of mutations)
barplot_sum4 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 10))+
  ggtitle(label = paste("Genetics for", classes, collapse = " ")) + theme(plot.title = element_text(color="black", size=10, face="bold"))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size = 10))

# png(filename = paste0("../results/9_genebarplot_",classes[1],".png"), width = 900, height = 500)
print(barplot_sum4)
# dev.off()


###5
classes <- c("ELN2022_intermediate")
mutation_select <- mutation_covariate_data[mutation_covariate_data$ELN2022_IPSSM %in% classes,]

subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL
# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
melted_data_sum$Cluster <- as.factor(melted_data_sum$Cluster)
levels(melted_data_sum$Cluster) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

# Stacked barplot (sum of mutations)
barplot_sum5 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 10))+
  ggtitle(label = paste("Genetics for", classes, collapse = " ")) + theme(plot.title = element_text(color="black", size=10, face="bold"))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size = 10))

#png(filename = paste0("../results/9_genebarplot_",classes[1],".png"), width = 900, height = 500)
print(barplot_sum5)
#dev.off()

###6 
classes <- c("ELN2022_favorable")
mutation_select <- mutation_covariate_data[mutation_covariate_data$ELN2022_IPSSM %in% classes,]

subdata <- mutation_select[,c(2:46,60)]
mutation_sum <- aggregate(subdata, by = list(subdata$`cluster_results$clustermembership`), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`cluster_results$clustermembership` <- NULL
# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
melted_data_sum$Cluster <- as.factor(melted_data_sum$Cluster)
levels(melted_data_sum$Cluster) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

# Stacked barplot (sum of mutations)
barplot_sum6 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 10))+
  ggtitle(label = paste("Genetics for", classes, collapse = " ")) + theme(plot.title = element_text(color="black", size=10, face="bold"))+
  theme(legend.text=element_text(size=8), legend.title = element_text(size = 10))

#png(filename = paste0("../results/9_genebarplot_",classes[1],".png"), width = 900, height = 500)
print(barplot_sum6)
#dev.off()


pdf(file = "./figures/SuppS6.pdf", width = 14, height = 17)
    ggarrange(barplot_sum6, barplot_sum5, barplot_sum4, barplot_sum1,barplot_sum2,barplot_sum3,ncol = 1,nrow = 6,labels = letters[1:6])
dev.off()



#### FIG SUPP 8

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
clusterAML <- readRDS("../comparison/AMLonly_6clus.rds") #readRDS("./AMLonly_9clus.rds")
clusterMDS <- readRDS("../comparison/MDSonly_4clus.rds") #clusterMDS <- readRDS("./MDSonly_6clus.rds")

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
levels(TOTmat$group) <-c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
combimembers <- merge(TOTmat,combimembers, by = "ID")

# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
totdata <- merge(survdata,combimembers, by="ID")

totdata$WHO_2016 <- as.factor(totdata$WHO_2016)
totdata$WHO_2022 <- as.factor(totdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(totdata$IPSSR_ELN)
totdata$ELN2022_IPSSM <- as.factor(survdata$ELN2022_IPSSM)
totdata$group.y <- as.factor(totdata$group.y)
#totdata$group.y <- factor(totdata$group.y, levels = c("1_AML", "2_AML", "3_AML", "4_AML", "5_AML", "1_MDS", "2_MDS", "3_MDS", "4_MDS", "5_MDS", "6_MDS"))

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
colours_clusters <- c("#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")

plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')

mycolor2 <- c("#202020","#771122","#AA4455","#DD7788","#CC99BB",
              "#114477","#4477AA","#77AADD", "#77CCCC","#88CCAA","#117744")


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")
library(viridis)
pos1 <- ggsurvplot(
  os,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = viridis(9), #mycolor2,#colours_clusters, # personalized colours
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
  risk.table.fontsize = 3,
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
  risk.table.fontsize = 3,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title="Overall Survival (n=7480)",
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)


pos2

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
  palette = viridis(9),#mycolor2[c(6:11)],#colours_clusters, # personalized colours
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
  risk.table.fontsize = 3,
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
  risk.table.fontsize = 3,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("Overall Survival - MARCS (n=", length(clinical$group), ")"),
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
  palette = plasma(9),#mycolor2[c(1:5)],#colours_clusters, # personalized colours
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
  risk.table.fontsize = 3,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("Overall Survival - AML only (n=", length(clinical$group), ")"),
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
  
)

amlos1

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
  risk.table.fontsize = 3,
  pval.size =4,
  pval.coord = c(0, 0.07),
  font.title= 16,
  title=paste0("Overall Survival - MARCS (n=", length(clinical$group), ")"),
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
  risk.table.fontsize = 3,
  pval.size =4,
  title=paste0("Time to AML transformation - MARCS (n=", length(clinical$group), ")"),
  font.title= 16,
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table 
  
)


##
clinical <- list()
#pfsdata <- pfsdata[grep(pattern = "MDS/AML", x = totdata$ICC),] # RESTRICT TO MDS/AML, comment for the whole dataset 
clinical$event <- pfsdata$AMLt_STATUS
clinical$time <- pfsdata$AMLt_YEARS
clinical$group <- pfsdata$group.y
clinical$groupSep <- pfsdata$group.y#pfsdata$ELN2022_IPSSM
clinical$groupTot <- pfsdata$group.x


# Analysis by Fritz on Jan 19th, 2024 / Update Marco 02.01.2024 - ANALYSIS OF PROGRESSION FREE SURVIVAL 

# add sex and age
clinical$sex <- pfsdata$SEX.x
clinical$age <- pfsdata$AGE.x

# add cancer type
clinical$is_aml <- rep(0,length(clinical$age))
clinical$is_aml[grep("AML",pfsdata$group.y)] <- 1


# Question One: How predictive are we beyond the clinical classifications (AML MDS classifications) and beyond the clinical data (age and sex)?

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  age + sex + groupSep + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + age + sex + groupSep + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(unique(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("TRANSFORMATION Combined vs Separated Total: ", paste(round(LR, 1), "&", pvalue, sep = " ")))

# Question Two: Which score is more predictive in a direct comparison? 

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(unique(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("Result for Combined Clustering: ", paste(round(LR, 1), "&", pvalue, sep = " ")))

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupSep + age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(unique(clinical$groupSep)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("Result for Separated Classification: ", paste(round(LR, 1), "&", pvalue, sep = " ")))


pfs <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")

pfs1 <- ggsurvplot(
  pfs,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = viridis(9),#mycolor2[c(6:11)], # personalized colours
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
  risk.table.fontsize = 3,
  pval.size =4,
  title=paste0("Time to AML transformation - MDS only (n=", length(clinical$group), ")"),
  font.title= 16,
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table 
  
)


# pdf(file = "../figures/Supp_compare_OStot.pdf",  width = 16.7, height = 8)
# ggarrange( pos2$plot, pos1$plot, pos2$table, pos1$table,
#            ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
# dev.off() 
# 
# pdf(file = "../figures/Supp_compare_OSAML.pdf",  width = 16.7, height = 8)
# ggarrange( amlos2$plot, amlos1$plot, amlos2$table, amlos1$table,
#            ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
# dev.off() 
# 
# pdf(file = "../figures/Supp_compare_Tot.pdf",  width = 16.7, height = 22)
# ggarrange( amlos2$plot, amlos1$plot, amlos2$table, amlos1$table,
#            mdsos2$plot, mdsos1$plot, mdsos2$table, mdsos1$table,
#            pfs2$plot, pfs1$plot, pfs2$table, pfs1$table,
#            ncol = 2, nrow = 6 , labels = c(letters[1:2], "", "",letters[3:4], "", "",letters[5:6],"",""), heights = c(1,0.4,1,0.4,1,0.4) )
# dev.off() 
# 
# pdf(file = "../figures/Supp_compare_OSMDS.pdf",  width = 16.7, height = 8)
# ggarrange( mdsos2$plot, mdsos1$plot, mdsos2$table, mdsos1$table,
#            ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
# dev.off() 
# 
# 
# pdf(file = "../figures/Supp_compare_PFS.pdf",  width = 16.7, height = 8)
# ggarrange( pfs2$plot, pfs1$plot, pfs2$table, pfs1$table,
#            ncol = 2, nrow = 2 , labels = c(letters[1:2], "", "",letters[3:4], "", ""), heights = c(1,0.4,1,0.4) )
# dev.off() 


pos1
pos2

amlos1
amlos2  

mdsos1
mdsos2

pfs1
pfs2



##### Analysis by Fritz on Jan 19th, 2024 / Update Marco 02.01.2024 - ANALYSIS OF OVERALL SURVIVAL 


clinical <- list()
#totdata <- totdata[grep(pattern = "MDS/AML", x = totdata$ICC),] # RESTRICT TO MDS/AML, comment for the whole dataset 
clinical$age <- totdata$AGE
clinical$sex <- totdata$SEX
clinical$event <- totdata$OS_STATUS
clinical$time <- totdata$OS
clinical$groupTot <- totdata$group.x
clinical$groupSep <- totdata$group.y#totdata$ELN2022_IPSSM
clinical$groupTot <- as.factor(clinical$groupTot)
clinical$groupSep <- as.factor(clinical$groupSep)

clinical$is_aml <- rep(0,length(clinical$age))
clinical$is_aml[grep("AML",totdata$group.y)] <- 1

# Question One: How predictive are we beyond the clinical classifications (AML MDS classifications) and beyond the clinical data (age and sex)?

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  age + sex + groupSep + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + age + sex + groupSep + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(unique(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("OVERALL SURVIVAL Combined vs Separated Total: ", paste(round(LR, 1), "&", pvalue, sep = " ")))

# Question Two: Which score is more predictive in a direct comparison? 

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(unique(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("Result for Combined Clustering: ", paste(round(LR, 1), "&", pvalue, sep = " ")))

tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]

groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupSep + age + sex + is_aml, data = clinical, na.action = "na.omit"))$logtest[1]
LR <- round((groupTest-tissueTest)/2, 1)
pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(unique(clinical$groupSep)) - 1, lower.tail = FALSE), 500)
# coxResults[1,1:2] <- c(LR,pvalue)
print(paste0("Result for Separated Classification: ", paste(round(LR, 1), "&", pvalue, sep = " ")))



# 
# 
# clinical <- list()
# clinical$age <- totaml$AGE
# clinical$sex <- totaml$SEX
# clinical$event <- totaml$OS_STATUS
# clinical$time <- totaml$OS
# clinical$groupTot <- totaml$group.x
# clinical$groupSep <- totaml$group.y
# clinical$groupTot <- as.factor(clinical$groupTot)
# clinical$groupSep <- as.factor(clinical$groupSep)
# 
# tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit"))$logtest[1]
# 
# groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + groupSep + age + sex, data = clinical, na.action = "na.omit"))$logtest[1]
# LR <- round((groupTest-tissueTest)/2, 1)
# pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# # coxResults[1,1:2] <- c(LR,pvalue)
# print(paste0("Combined vs Separated AML: ", paste(round(LR, 1), "&", pvalue, sep = " ")))



# #### CONCORDANCE INDEX
# 
# library(survival)
# fitSep <- coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit")
# fitTot <- coxph(Surv(time, as.numeric(event)) ~  groupTot+ age + sex, data = clinical, na.action = "na.omit")
# 
# concordance(fitSep)
# concordance(fitTot)
# 
# clinical <- list()
# clinical$age <- totmds$AGE
# clinical$sex <- totmds$SEX
# clinical$event <- totmds$OS_STATUS
# clinical$time <- totmds$OS
# clinical$groupTot <- totmds$group.x
# clinical$groupSep <- totmds$group.y
# clinical$groupTot <- as.factor(clinical$groupTot)
# clinical$groupSep <- as.factor(clinical$groupSep)
# 
# tissueTest <- summary(coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit"))$logtest[1]
# 
# groupTest <- summary(coxph(Surv(time, as.numeric(event)) ~ groupTot + groupSep + age + sex, data = clinical, na.action = "na.omit"))$logtest[1]
# LR <- round((groupTest-tissueTest)/2, 1)
# pvalue <- round(pchisq(q = groupTest-tissueTest, df = length(levels(clinical$groupTot)) - 1, lower.tail = FALSE), 500)
# # coxResults[1,1:2] <- c(LR,pvalue)
# print(paste0("Combined vs Separated MDS: ", paste(round(LR, 1), "&", pvalue, sep = " ")))
# 
# library(survival)
# fitSep <- coxph(Surv(time, as.numeric(event)) ~  groupSep+ age + sex, data = clinical, na.action = "na.omit")
# fitTot <- coxph(Surv(time, as.numeric(event)) ~  groupTot+ age + sex, data = clinical, na.action = "na.omit")
# 
# print("Separated MDS: ")
# concordance(fitSep)
# print("Combined MDS: ")
# concordance(fitTot)
# 

### 
pdf("../analysis/figures/SupplS8.pdf", width = 20, height = 25)
ggarrange(pos2$plot,pos1$plot,pos2$table,pos1$table,amlos2$plot,amlos1$plot,amlos2$table,amlos1$table,mdsos2$plot,mdsos1$plot,mdsos2$table,mdsos1$table,pfs2$plot,pfs1$plot,pfs2$table,pfs1$table, 
          ncol = 2, nrow = 8, heights = c(2,1), labels = c("a", "b","", "", "c", "d", "", "", "e",  "f", "", "", "g", "h", "", ""))
dev.off()


##### SUPPL 13 

library(clustNet)
library(parallel)
library(xlsx)

rm(list=ls())

cluster_results <- readRDS("../results/euler_memberships_8k_9clusters.rds")


validata <- read.table("../validation/mda_binary-mutationCovariate-matrix.txt")[,-58]
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
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
levels(melted_data_sum$Cluster) <- lab # LETTERS[1:length(unique(cluster_results$clustermembership))]
# 
melted_data_sum$Cluster[melted_data_sum$Cluster=="1"] <- "UHR"
melted_data_sum$Cluster[melted_data_sum$Cluster=="2"] <- "HR1"
melted_data_sum$Cluster[melted_data_sum$Cluster=="3"] <- "NPM1"
melted_data_sum$Cluster[melted_data_sum$Cluster=="4"] <- "HR2"
melted_data_sum$Cluster[melted_data_sum$Cluster=="5"] <- "INT2"
melted_data_sum$Cluster[melted_data_sum$Cluster=="6"] <- "HR3"
melted_data_sum$Cluster[melted_data_sum$Cluster=="7"] <- "INT1"
melted_data_sum$Cluster[melted_data_sum$Cluster=="8"] <- "LR1"
melted_data_sum$Cluster[melted_data_sum$Cluster=="9"] <- "LR2"
# 

colours_clusters <- c("#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")
# Stacked barplot (sum of mutations)
barplot_sum0 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
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
levels(clinical$group) <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

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

colours_clusters <- c("#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")

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
  #legend.labs = sort(as.character(unique(filtered_clinical_aml$group))), 
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
  #legend.labs = sort(as.character(unique(filtered_clinical_mds$group))), 
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


pdf(file = "./figures/SupplS13.pdf", width = 16.7, height = 20)
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

##### SUPPLEMENTARY S9

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

cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")
#clusterAML <- readRDS("./MDSonly_6clus.rds")
#clusterMDS <- readRDS("./MDSonly_6clus.rds")
#readRDS("../results/euler_memberships_8k_9clusters.rds")

# remove this line in the future (after package update on euler)
# cluster_res$clustermembership <- cluster_res$clustermembership[[1]]
# cluster_res <- readRDS("../results/custer_res.rds")

table(cluster_res$clustermembership)

# mut_cov_data <- read.csv("../data/undivided_binary_matrix.csv")

# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
survdata <- survdata[survdata$ID %in% rownames(cluster_res$data),]
clumatrix <- data.frame(ID=rownames(cluster_res$data), group=cluster_res$clustermembership)
survdata <- merge(survdata,clumatrix, by = "ID")

# survdata$last_fu <- as.Date(survdata$last_fu, format = "%d.%m.%Y")
# survdata$firstdiagnosis_date <- as.Date(survdata$firstdiagnosis_date, format = "%d.%m.%Y")
# survdata$relapse_date <- as.Date(survdata$relapse_date, format = "%d.%m.%Y")
# survdata$type <- as.factor(survdata$)
survdata$WHO_2016 <- as.factor(survdata$WHO_2016)
survdata$WHO_2022 <- as.factor(survdata$WHO_2022)
#survdata$IPSSR_ELN <- as.factor(survdata$IPSSR_ELN)
survdata$ELN2022_IPSSM <- as.factor(survdata$ELN2022_IPSSM)
survdata$Cluster <- survdata$group
survdata$Cluster <- as.factor(survdata$Cluster)
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(survdata$Cluster) <- lab
survdata <- survdata[grep("MDS/AML",survdata$ICC),]

survdata <- survdata[!survdata$ELN2022_IPSSM=="IPSSM_NA",]

ipssm_vector <- c("IPSSM_Low", "IPSSM_Moderate-Low", "IPSSM_High", "IPSSM_Very-Low", "IPSSM_Moderate-High", "IPSSM_Very-High")


p <- c()

for(i in ipssm_vector){
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
  colours_clusters <- c("blue", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")#"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")
  
  #plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
  #legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')
  mypal=colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])]
  dimensione=length(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])])-1
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
  library("survminer")
  file<-paste0("../comparison/Supp_","MDSAML_",i,".png")
  png(filename = file, width = 2300, height = 2000, res = 280)
  p[i] <- ggsurvplot(
    os,                     # survfit object with calculated statistics.
    data = clinical,             # data used to fit survival curves.
    #palette = colours_clusters,#colours_clusters, # personalized colours
    palette = alpha(mypal,c(1,rep(0.55, dimensione))), 
    # legend.labs = sort(as.character(unique(survdata$Cluster))),
    legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
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
    pval.coord = c(0, 0.07),
    title= paste0("OS for MDS/AML - ", gsub(pattern = "_", replacement = " ", i), " (n=",length(clinical$event)/2,")"),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
    
  )
  print(p[i])
  
  dev.off()
}

# pdf(file = "../comparison/Supplementary_MDSAML_IPSSM_OS.pdf", width = 11, height = 10)
# ggarrange(p$IPSSM_Low, ncol = 2, nrow = 2 , labels = letters[1:4])
# dev.off() 

#

pdf(file = "./figures/SuppS9.pdf", width = 16.7, height = 10)
ggarrange(p$`IPSSM_Very-High`,p$IPSSM_High,
          p$`IPSSM_Moderate-High`,p$`IPSSM_Moderate-Low`, 
          p$IPSSM_Low,p$`IPSSM_Very-Low`, ncol = 3, nrow = 2 , labels = letters[1:6])
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
#cluster_res <- cluster_res <- clusterAML <- readRDS("./MDSonly_6clus.rds")
#clusterMDS <- readRDS("./MDSonly_6clus.rds")

# remove this line in the future (after package update on euler)
# cluster_res$clustermembership <- cluster_res$clustermembership[[1]]
# cluster_res <- readRDS("../results/custer_res.rds")

table(cluster_res$clustermembership)

# mut_cov_data <- read.csv("../data/undivided_binary_matrix.csv")

# import survival data
survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
survdata <- survdata[survdata$ID %in% rownames(cluster_res$data),]
clumatrix <- data.frame(ID=rownames(cluster_res$data), group=cluster_res$clustermembership)
survdata <- merge(survdata,clumatrix, by = "ID")

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
survdata <- survdata[grep("MDS/AML",survdata$ICC),]
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(survdata$Cluster) <- lab #LETTERS[1:length(unique(survdata$Cluster))]

survdata <- survdata[!survdata$ELN2022_IPSSM=="IPSSM_NA",]

p <- c()

#,"IPSSM_Low","IPSSM_Moderate-Low"

for(i in c("IPSSM_High","IPSSM_Very-High","IPSSM_Moderate-High")){
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
  colours_clusters <- c("blue", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73") #"#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")
  
  #plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
  #legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')
  mypal=colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])]
  dimensione=length(colours_clusters[c(1, c(seq(1:9)+1)[as.character(sort(unique(survdata$Cluster))) %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])])-1
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
  library("survminer")
  file<-paste0("../comparison/Supp_","MDSAML_",i,"_AMLt.png")
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
    ylab = "Transformation probability",
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.fontsize = 3,
    pval.size =4,
    pval.coord = c(0, 0.07),
    title= paste0("AMLt for MDS/AML - ", gsub(pattern = "_", replacement = " ", i), " (n=",length(clinical$event)/2,")"),
    risk.table.y.text = FALSE, # show bars instead of names in text annotations
    legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)]))# in legend of risk table
    
  )
  print(p[i])
  
  dev.off()
}
pdf(file = "./figures/SuppS10.pdf", width = 16.7, height = 5)
ggarrange(p$`IPSSM_Very-High`,p$IPSSM_High,
          p$`IPSSM_Moderate-High`, ncol = 3, nrow = 1 , labels = letters[1:4])
dev.off() 



#### SUPPLEMENTARY S19 - now S15

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


# clinical$age <- #here needs to go the age of the patients#



### ALL
groups <- data.frame(cluster_res$clustermembership)
groups$ID <- rownames(cluster_res$data)
groups$cluster_res.clustermembership <- as.factor(groups$cluster_res.clustermembership)
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(groups$cluster_res.clustermembership) <- lab#LETTERS[1:length(unique(clinical$group))]
subdata53 <-  merge(survdata, groups, by="ID")

quartos <- c()
for(i in lab){
  savedata <- subdata53
  subdata53 <- subdata53[subdata53$cluster_res.clustermembership==i,]
  subdata53 <- subdata53 %>% mutate(blast_quartile = ntile(BM_BLASTS, 4))
  subdata53 <- subdata53 %>% mutate(median = ntile(BM_BLASTS, 2))
  #subdata53 <- subdata53[subdata53$blast_quartile %in% c(1,4),]
  clinical <- list()
  clinical$event <- subdata53$OS_STATUS
  clinical$time <- subdata53$OS
  clinical$quartile <- subdata53$blast_quartile
  clinical$median <- subdata53$median
  clinical$type <- subdata53$WHO_2016
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ quartile, data = clinical)
  library("survminer")
  #png("../results/9_TP53_mut_KM.png", width = 1200, height = 1000, res=150)
  quartos[i] <- ggsurvplot(
    os,                     # survfit object with calculated statistics.
    data = clinical,             # data used to fit survival curves.
    #palette = colours_clusters, # personalized colours
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = F,         # show confidence intervals for 
    # point estimates of survival curves.
    xlim = c(0,10),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in years",   # customize X axis label.
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.fontsize = 3,
    pval.size =4,
    title=paste0("OS in group ",i, " per blast quartiles"),
    legend.labs=c("1st quart.", "2nd quart.","3rd quart.","4th quart."),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )  
  asd <- coxph(Surv(OS, OS_STATUS) ~ BM_BLASTS, data = subdata53)
  print(i)
  print(asd)
  subdata53 <- savedata
}


# Assuming your data frame is named df and the column of interest is named condition
# clinical <- as.data.frame(clinical)
# clinical <- clinical %>%
#   mutate(type = ifelse(grepl("AML", type), "AML", "MDS"))
# os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ quartile, data = clinical)

typos <- c()
for(i in lab){
  savedata <- subdata53
  subdata53 <- subdata53[subdata53$cluster_res.clustermembership==i,]
  subdata53 <- subdata53 %>% mutate(blast_quartile = ntile(BM_BLASTS, 4))
  subdata53 <- subdata53 %>% mutate(median = ntile(BM_BLASTS, 2))
  #subdata53 <- subdata53[subdata53$blast_quartile %in% c(1,4),]
  clinical <- list()
  clinical$event <- subdata53$OS_STATUS
  clinical$time <- subdata53$OS
  clinical$quartile <- subdata53$blast_quartile
  clinical$median <- subdata53$median
  clinical$type <- subdata53$WHO_2016
  clinical <- as.data.frame(clinical)
  clinical <- clinical %>%  mutate(type = ifelse(grepl("AML", type), "AML", "MDS"))
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ clinical$type, data = clinical)
  typos[i] <- ggsurvplot(
    os,                     # survfit object with calculated statistics.
    data = clinical,             # data used to fit survival curves.
    #palette = colours_clusters, # personalized colours
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    conf.int = F,         # show confidence intervals for 
    # point estimates of survival curves.
    xlim = c(0,10),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in years",   # customize X axis label.
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.fontsize = 3,
    pval.size =4,
    title=paste0("OS in group ",i, " per disease type"),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  ); p
  print(i)
  print("AML - 5yrs OS")
  print(sum(as.integer(clinical$time[clinical$type=="AML"] >5))/length(clinical$time[clinical$type=="AML"]))

  print("MDS - 5yrs OS")
  print(sum(as.integer(clinical$time[clinical$type=="MDS"] >5))/length(clinical$time[clinical$type=="MDS"]))

  subdata53 <- savedata
}


png("./figures/SuppS15.png", width = 3600, height = 4500, res = 250)
ggarrange(typos$UHR,typos$HR1,typos$HR2,typos$HR3, typos$NPM1, typos$INT1, typos$INT2, typos$LR1, typos$LR2,
          heights =  c(1,1), labels = c("a","b","c", "d","e", "f","g","h","i"), ncol = 3, nrow = 3)
dev.off()

png("./figures/SuppS16.png", width = 3600, height = 4500, res = 250)
ggarrange(quartos$UHR,quartos$HR1,quartos$HR2,quartos$HR3, quartos$NPM1, quartos$INT1, quartos$INT2, quartos$LR1, quartos$LR2,
          heights =  c(1,1), labels = c("a","b","c", "d","e", "f","g","h","i"), ncol = 3, nrow = 3)
dev.off()





######SUPPLEMENTARY S9 - tAML in MDS/AMl 

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
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
levels(survdata$Cluster) <- lab
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
  # 
  colours_clusters <- c("blue", "#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")# "#000000","#772266","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")
  
  #plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
  #legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')
  mypal=colours_clusters[c(1, c(seq(1:9)+1)[lab %in% unique(clinical$group)[grep(x = unique(clinical$group), pattern = i, invert = T)]])]
  dimensione=length(mypal)-1
  
  
  
  os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
  library("survminer")
  #file<-paste0("../figures/Supp_",gsub("MDS/AML ", "MDSAML_",i),".png")
  #png(filename = file, width = 2300, height = 2000, res = 280)
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
    legend.labs = c(unique(clinical$group)[1], sort(unique(clinical$group)[2:length(clinical$group)])),
    title= paste0("Time to AML transformation - ", gsub(pattern = "_", replacement = " ", i), " (n=",length(clinical$event)/2,")"),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
    
  )
  print(p[i])
  
  #dev.off()
}

pdf(file = "./figures//SuppS9_tAML.pdf", width = 16.7, height = 10)
ggarrange(p$`MDS/AML MR-Cyto`,p$`MDS/AML MR-Gen`, p$`MDS/AML TP53`, p$`MDS/AML NOS`, ncol = 2, nrow = 2 , labels = letters[1:4])
dev.off() 

##### ALLUVIAL 

library(ggalluvial)
library(graphClust)

# Prepare session, load packages
rm(list=ls())
library(survival)
library(RColorBrewer)

# Load classification by mutation profile (cluster assignment) # CLUSTERS 
cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")

table(cluster_res$clustermembership)
mutation_covariate_data <- read.csv("../data/undivided_binary_only_matrix_8k.csv")

survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
survdata$clustermembership <- cluster_res$clustermembership
survdata$clustermembership <- as.factor(survdata$clustermembership)
levels(survdata$clustermembership) <- LETTERS[1:length(unique(cluster_res$clustermembership))]

survdata <- survdata[!survdata$ELN2022_IPSSM == "IPSSM_NA",]
survdata <- survdata[!survdata$IPSSR_ELN == "IPSSR_NA",]

# create alluvial plot
# 
# allu <- ggplot(data = survdata,
#                aes(axis2 = WHO_2016, axis3 = WHO_2022, axis4 = ICC, axis1 = ELN2022_IPSSM )) +
#   geom_alluvium(aes(fill = ELN2022_IPSSM)) +
#   geom_stratum(width=0.50) +
#   geom_text(stat = "stratum", 
#             aes(label = after_stat(stratum)),size= 3) +
#   scale_x_discrete(limits = c("WHO_2016", "WHO_2022", "ICC", "ELN2022_IPSSM"),
#                    expand = c(0.1, 0.15, 0.1, 0.05)) +
#   # theme_void() +
#     theme(legend.position = "none") ; allu


allu <- ggplot(data = survdata,
               aes(axis2 = WHO_2016, axis3 = WHO_2022, axis4 = ICC, axis1 = ELN2022_IPSSM)) +
  geom_alluvium(aes(fill = ELN2022_IPSSM)) +
  geom_stratum(width = 0.50) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("ELN2022_IPSSM","WHO_2016", "WHO_2022", "ICC"),
                   expand = c(0.1, 0.15, 0.1, 0.05)) +
  theme_void() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 13, vjust = 2), # Adjust the margin to reduce space
    plot.margin = unit(c(-1, -.5, 0, -.5), "lines")); allu # Customize the x-axis text appearance


# allu <- ggplot(data = survdata,
#                aes(axis2 = WHO_2016, axis3 = WHO_2022, axis4 = ICC, axis1 = ELN2022_IPSSM)) +
#   geom_alluvium(aes(fill = ELN2022_IPSSM)) +
#   geom_stratum(width = 0.50) +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3, 
#              position = position_nudge(x = -0.25), # Nudge labels to the side
#              hjust = 0, vjust = 0.5, # Adjust horizontal and vertical justification
#              label.size = NA, # Remove label border
#              fill = "white") + # White background for labels
#   scale_x_discrete(limits = c("ELN2022_IPSSM","WHO_2016", "WHO_2022", "ICC"),
#                    expand = c(0.1, 0.15, 0.1, 0.05)) +
#   theme_void() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(angle = 90, hjust = 0.6))
# allu


png(filename = "./figures/Supp_Alluvial.png", width = 2400, height = 3000, res = 150)
allu
dev.off()

