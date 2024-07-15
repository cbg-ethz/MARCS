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
colours_clusters2 <- c("#2f0000","#BB3E03","#CA6702","#EE9B00","#0A9396","#005F73") #c("#000000","#772266ff","#117777","#771122aa","#114477","#CC99BBcc","#88CCAA","#117744","#77AADD")


plot(survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical), col=colourysdots, ylab='Survival probability', xlab='Survival (years)', mark.time = T)
legend(x = 12.9, y = 1.0, legend = paste(levels(clinical$group), sep=' '), pch = 15, col=colourysdots, cex=1, ncol = 2,title="Cluster", bg='white')

mycolor2 <- c("#202020","#771122","#AA4455","#CC99BB","#774411","#AA7744",
              "#DDAA77","#DD7788","#DDCC77","#77CCCC","#114477","#4477AA",
              "#77AADD","#AA4488","#CC99BB","#DDCC77","#777711","#AAAA44","#117744","#44AA77",
              "#88CCAA","#1122AA","#44AA00","#77CC77", "#AAAA44","#DDDD77","blue", "red", "yellow") ## COLOURS NEEDS TO BE ADJUSTED 


os <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")


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
lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
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

pfs <- survfit(Surv(time = as.numeric(time), event = as.numeric(event)) ~ group, data = clinical)
library("survminer")


oslist <- c()
for(i in 1:length(colours_clusters)){
  asd <- c("gray90","gray91","gray92","gray93","gray94","gray95","gray96","gray97","gray98")
  asd[i] <- colours_clusters[i]
  
  
  oslist[lab[i]] <- ggsurvplot(
    os,                     # survfit object with calculated statistics.
    data = clinical,             # data used to fit survival curves.
    palette = asd, #colours_clusters, # personalized colours
    risk.table = F,       # show risk table.
    pval = F,             # show p-value of log-rank test.
    conf.int = T,         # show confidence intervals for 
    # point estimates of survival curves.
    xlim = c(0,10),         # present narrower X axis, but not affect
    # survival estimates.
    xlab = "Time in years",   # customize X axis label.
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_light(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.fontsize = 4,
    pval.size =4,
    pval.coord = c(0, 0.15),
    font.title= 24,
    title="Overall Survival",
    legend.labs = lab, # in legend of risk table,
    font.legend = list(size = 16),
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )

  ggpar(oslist[lab[i]], 
        font.main = c(24, "bold"),
        font.x = c(14, "bold"),
        font.y = c(14, "bold"),
        font.caption = c(14, "bold"), 
        font.legend = c(16, "bold"), 
        font.tickslab = c(14, "bold"))

}



lab2 <- c("UHR", "HR2", "HR3", "INT1", "LR1", "LR2")
pfslist <- c()
for(i in 1:length(colours_clusters2)){
  asd2 <- c("gray90","gray91","gray92","gray93","gray94","gray95")
  asd2[i] <- colours_clusters2[i]
  
  
  pfslist[lab2[i]] <- ggsurvplot(
      pfs,                     # survfit object with calculated statistics.
      data = clinical,             # data used to fit survival curves.
      palette = asd2, #colours_clusters, # personalized colours
      risk.table = F,       # show risk table.
      pval = F,             # show p-value of log-rank test.
      conf.int = T,         # show confidence intervals for 
      # point estimates of survival curves.
      xlim = c(0,10),         # present narrower X axis, but not affect
      # survival estimates.
      ylab = "Transformation probabiliy",
      xlab = "Time in years",   # customize X axis label.
      break.time.by = 1,     # break X axis in time intervals by 500.
      ggtheme = theme_light(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.fontsize = 4,
      pval.size =4,
      title="Time to AML transformation",
      font.title= 24,
      font.legend = list(size = 16),
      legend.labs = c("UHR","HR2","HR3", "INT1","LR1","LR2"), # in legend of risk table,
      risk.table.y.text = FALSE # show bars instead of names in text annotations
      # in legend of risk table 
    )
}


asd2 <- c("gray90","gray91","gray92","gray93","gray94","gray95")
pfsblank <- ggsurvplot(
  pfs,                     # survfit object with calculated statistics.
  data = clinical,             # data used to fit survival curves.
  palette = asd2, #colours_clusters, # personalized colours
  risk.table = F,       # show risk table.
  pval = F,             # show p-value of log-rank test.
  conf.int = T,         # show confidence intervals for 
  # point estimates of survival curves.
  xlim = c(0,10),         # present narrower X axis, but not affect
  # survival estimates.
  ylab = "Transformation probabiliy",
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = F, # colour risk table text annotations.
  risk.table.fontsize = 4,
  pval.size =4,
  title="Time to AML transformation",
  font.title= 24,
  font.legend = list(size = 16),
  legend.labs = c("UHR","HR2","HR3", "INT1","LR1","LR2"), # in legend of risk table,
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table 
)


# ggpar(pfsblank, 
#       font.main = c(24, "bold"),
#       font.x = c(14, "bold"),
#       font.y = c(14, "bold"),
#       font.caption = c(14, "bold"), 
#       font.legend = c(16, "bold"), 
#       font.tickslab = c(14, "bold"))


png("../shiny/www/kaplan_maier/KM_UHR.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$UHR, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfslist$UHR, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
           ncol = 2, nrow = 1 )
dev.off()

png("../shiny/www/kaplan_maier/KM_HR1.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$HR1, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfsblank$plot, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()



png("../shiny/www/kaplan_maier/KM_HR2.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$HR2, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfslist$HR2, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()


png("../shiny/www/kaplan_maier/KM_HR3.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$HR3, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfslist$HR3, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()


png("../shiny/www/kaplan_maier/KM_INT1.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$INT1, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfslist$INT1, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()


png("../shiny/www/kaplan_maier/KM_INT2.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$INT2, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfsblank$plot, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()

png("../shiny/www/kaplan_maier/KM_NPM1.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$NPM1, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfsblank$plot, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()


png("../shiny/www/kaplan_maier/KM_LR1.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$LR1, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfslist$LR1, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()


png("../shiny/www/kaplan_maier/KM_LR2.png", width = 1920, height = 800)
ggarrange(
  ggpar(oslist$LR2, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ggpar(pfslist$LR2, 
        font.main = c(24, "bold"),
        font.x = c(20),
        font.y = c(20),
        font.caption = c(16, "bold"), 
        font.legend = c(20, "bold"), 
        font.tickslab = c(18, "bold")),
  ncol = 2, nrow = 1 )
dev.off()
