library(ggalluvial)
library(graphClust)

# Prepare session, load packages

library(survival)
library(RColorBrewer)

# get bar plots of cluster assignment

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
mutation_covariate_data <- mutation_covariate_data[!mutation_covariate_data$ELN2022_IPSSM=="IPSSM_NA",] 


### BARPLOT ON OLD CANCER STRATIFICATION

k_clust <- length(levels(as.factor(cluster_results$clustermembership)))

n_cancer_type <- length(table(mutation_covariate_data$IPSSR_ELN))

cancer_table <- matrix(NA, k_clust, n_cancer_type)
cancer_names <- names(table(mutation_covariate_data$IPSSR_ELN))

colnames(cancer_table) <- cancer_names
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
rownames(cancer_table) <- lab #LETTERS[1:k_clust]


for (ii in 1:k_clust){
  temp_count <- c()
  for (jj in 1:n_cancer_type){
    temp_count <- c(temp_count, length(which(mutation_covariate_data$IPSSR_ELN[which(cluster_results$clustermembership==ii)]==cancer_names[jj])))
  }
  cancer_table[ii,] <- temp_count
}

mycolor <- c("#DD7788", "#771122", "#DDDD77", "#117777", "#774411","#AA7744",
             "#DDAA77","#777711","#AAAA44","#DDDD77","#117744","#44AA77",
             "#88CCAA","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
             "#77AADD","#771155","#AA4488","#CC99BB")
mycolor <- c("#202020","#771122","#AA4455","#DD7788","#774411","#AA7744",
             "#DDAA77","#117777","#44AAAA","#77CCCC","#114477","#4477AA",
             "#77AADD","#771155","#AA4488","#CC99BB")

mycolor2 <- c("#202020","#AA7744","#774411","#AA4455","#CC99BB","#ee7688dd",
              "#771122","#DDCC77","#77CCCC","#114477","#4477AA",
              "#77AADD","#AA4488","red3","#DDCC77","#777711","#AAAA44","#117744","#44AA77",
              "#88CCAA","#1122AA","#44AA00","#77CC77", "#771155","#AA4488","darkgrey")



### BARPLOT 

k_clust <- length(levels(as.factor(cluster_results$clustermembership)))

n_cancer_type <- length(table(mutation_covariate_data$ELN2022_IPSSM))

cancer_table <- matrix(NA, k_clust, n_cancer_type)
cancer_names <- names(table(mutation_covariate_data$ELN2022_IPSSM))

colnames(cancer_table) <- cancer_names
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
rownames(cancer_table) <- lab #LETTERS[1:k_clust]


for (ii in 1:k_clust){
  temp_count <- c()
  for (jj in 1:n_cancer_type){
    temp_count <- c(temp_count, length(which(mutation_covariate_data$ELN2022_IPSSM[which(cluster_results$clustermembership==ii)]==cancer_names[jj])))
  }
  cancer_table[ii,] <- temp_count
}


mycolor2 <- c("#2a204f","#AA7744","#774411","#AA4455","#CC99BB","#ee7688dd",
              "#771122","#DDCC77","#77CCCC","#114477","#4477AA",
              "#77AADD","#AA4488","red3","#DDCC77","#777711","#AAAA44","#202020","#44AA77",
              "#88CCAA","#1122AA","#44AA00","#77CC77", "#771155","#AA4488","darkgrey")


riskcolor <- c("#202020","#774411","#DDAA77","#440154", "#453581","#34618D","#24878E","#88CCAA","#77CC77")

# p1_riskNew <- ggbarplot(melt(cancer_table), "Var1", "value",
#                         fill = "Var2", palette = riskcolor, color = NA,
#                         label = TRUE, lab.col = NA)+ 
#   # label = TRUE, lab.col = "white", lab.pos = "in")+ 
#   scale_fill_manual(breaks=c("ELN2022_adverse","ELN2022_intermediate","ELN2022_favorable","IPSSM_Very-High", "IPSSM_High", "IPSSM_Moderate-High","IPSSM_Moderate-Low","IPSSM_Low", "IPSSM_Very-Low"), values=riskcolor)+
#   #scale_fill_manual(values=c("ELN2022_adverse"="#202020","ELN2022_intermediate"="#774411","ELN2022_favorable"="#DDAA77","IPSSM_Very-High"="#440154FF", "IPSSM_High"="#453581FF", "IPSSM_Moderate-High"="#34618DFF","IPSSM_Moderate-Low"="#24878EFF","IPSSM_Low" ="#40BC72FF", "IPSSM_Very-Low"="#CBE11EFF"))+
#   xlab("Clusters") + 
#   ylab("Number of Patients") +
#   theme(legend.key.size = unit(0.2, 'cm'),
#         legend.text = element_text(size=8)) +
#   
#   guides(fill=guide_legend(title="Risk Stratification (ELN2022,IPSS-M)"),col = FALSE)+
#   theme(legend.position="right"); p1_riskNew
# 
# p1_riskNew


### BARPLOT ON WHO2016 DIAGNOSIS 

k_clust <- length(levels(as.factor(cluster_results$clustermembership)))

n_cancer_type <- length(table(mutation_covariate_data$WHO_2016))

cancer_table <- matrix(NA, k_clust, n_cancer_type)
cancer_names <- names(table(mutation_covariate_data$WHO_2016))

colnames(cancer_table) <- cancer_names
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
rownames(cancer_table) <- lab #LETTERS[1:k_clust]

for (ii in 1:k_clust){
  temp_count <- c()
  for (jj in 1:n_cancer_type){
    temp_count <- c(temp_count, length(which(mutation_covariate_data$WHO_2016[which(cluster_results$clustermembership==ii)]==cancer_names[jj])))
  }
  cancer_table[ii,] <- temp_count
}


col2016 <- c("grey80","#ee4455ff","#CC99BB","#ee7688ff","#CC4477",
             "#771122","#771155","#aa0066",
             "#a03020", "#BB4455ff", "#991122","#771188", 
             "grey35",
             "#DDCC77","#777711",
             "#88CCAA","#117777","#44AAAA","#77CCCC",
             "#77AADD","#114477","#4477AA","#1177AA","darkblue",
             "#AAAA4477","#E4795EFF","#EA9B3D","darkgrey")




p1_2016 <- ggbarplot(melt(cancer_table), "Var1", "value",
                     fill = "Var2", color = NA, palette = col2016,
                     label = TRUE, lab.col = NA)+ 
  scale_fill_manual(breaks=c("aCML", "AML CEBPAbi","AML inv(16)", "AML NPM1","AML t(8;21)",
                             "AML inv(3)","AML MRC", "AML NOS", 
                             "AML RUNX1","AML t(15;17)", "AML t(6;9)", "AML t(9;11)",    
                             "CMML",
                             "HR-MDS", "LR-MDS", 
                             "MDS-del5q", "MDS-RS-MLD", "MDS-RS-SLD", "MDS-RS-SLD/MLD", 
                             "MDS-SLD","MDS-SLD/MLD","MDS-MLD", "MDS-EB1", "MDS-EB2",
                             "MDS-U","MDS/MPN-RS-T","MDS/MPN-U","other") , values = col2016)+
  # label = TRUE, lab.col = "white", lab.pos = "in")+ 
  xlab("Clusters") + 
  ylab("Number of Patients") +
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=12)) +
  # ylab("") +
  guides(fill=guide_legend(title="Type (WHO2016)", ncol=2),col = FALSE)+
  theme(legend.position="right")+
  theme(legend.text = element_text(size=12))

p1_2016


### BARPLOT ON WHO2022 DIAGNOSIS 

k_clust <- length(levels(as.factor(cluster_results$clustermembership)))

n_cancer_type <- length(table(mutation_covariate_data$WHO_2022))

cancer_table <- matrix(NA, k_clust, n_cancer_type)
cancer_names <- names(table(mutation_covariate_data$WHO_2022))

colnames(cancer_table) <- cancer_names
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
rownames(cancer_table) <- lab #LETTERS[1:k_clust]


for (ii in 1:k_clust){
  temp_count <- c()
  for (jj in 1:n_cancer_type){
    temp_count <- c(temp_count, length(which(mutation_covariate_data$WHO_2022[which(cluster_results$clustermembership==ii)]==cancer_names[jj])))
  }
  cancer_table[ii,] <- temp_count
}


#mycolor2 <- c("#2a204f","#AA7744","#774411","#AA4455","#CC99BB","#ee7688dd",
#               "#771122","#DDCC77","#77CCCC","#114477","#4477AA",
#               "#77AADD","#AA4488","red3","#DDCC77","#777711","#AAAA44","lightblue","#44AA77",
#               "#88CCAA","#1122AA","#44AA00","#77CC77", "#771155","#AA4488","darkgrey")

library(viridis)
#mycolor2 <- c(plasma(11),"#202020",viridis(6))
mycolor2 <- c("grey80","#ee4455ff","#CC99BB", "#771122","#771155","#aa0066","#ee7688ff","#CC4477","#BB4455ff","#a03020",  "#CC4477","grey35","#88CCAA","#1177AA","darkblue","#117777","#119944","#E4795EFF","#EA9B3D","#77AADD","#117744")

p1_2022 <- ggbarplot(melt(cancer_table), "Var1", "value",
                     fill = "Var2", color = "Var2", palette = mycolor2,
                     label = TRUE, lab.col = NA)+ 
  # label = TRUE, lab.col = "white", lab.pos = "in")+ 
  xlab("Clusters") + 
  ylab("Number of Patients") +
  # ylab("") +
  theme(legend.key.size = unit(0.4, 'cm'),
    legend.text = element_text(size=12)) +
  guides(fill=guide_legend(title="Type (WHO2022)", ncol=2),col = FALSE)+
  theme(legend.position="right")+
  theme(legend.text = element_text(size=12))

p1_2022



### BARPLOT ICC

k_clust <- length(levels(as.factor(cluster_results$clustermembership)))

n_cancer_type <- length(table(mutation_covariate_data$ICC))

cancer_table <- matrix(NA, k_clust, n_cancer_type)
cancer_names <- names(table(mutation_covariate_data$ICC))

colnames(cancer_table) <- cancer_names
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
rownames(cancer_table) <- lab #LETTERS[1:k_clust]


for (ii in 1:k_clust){
  temp_count <- c()
  for (jj in 1:n_cancer_type){
    temp_count <- c(temp_count, length(which(mutation_covariate_data$ICC[which(cluster_results$clustermembership==ii)]==cancer_names[jj])))
  }
  cancer_table[ii,] <- temp_count
}

# mycolor2 <- c("#2a204f","#AA7744","#774411","#AA4455","#CC99BB","#ee7688dd",
#               "#771122","#DDCC77","#77CCCC","#114477","#4477AA",
#               "#77AADD","#AA4488","lightblue","#DDCC77","#777711","#AAAA44","red3","#44AA77",
#               "#88CCAA","#1122AA","#44AA00","#77CC77", "#771155","#AA4488","darkgrey")

#mycolor2 <- c(plasma(9),"#202020",viridis(9))

mycolor2 <- c("grey80","#ee4477ff","#CC99BB", "#771122","#7500c0","#9900aa","#aa0066","#ee7688ff","#aa0066","#CC4477","#BB4455ff","#dd89aa", "#bb3020","#100000","grey35", "#88CCAA","darkblue","#117777","#119944","#aa0077","#770050","#500077","#006608","#550001","#E4795EFF","#EA9B3D","red","blue", "yellow", "green", "pink","darkred","grey")

p1_icc <- ggbarplot(melt(cancer_table), "Var1", "value",
                    fill = "Var2", color = "Var2", palette = mycolor2,
                    label = TRUE, lab.col = NA)+ 
  # label = TRUE, lab.col = "white", lab.pos = "in")+ 
  xlab("Clusters") + 
  ylab("Number of Patients") +
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=12)) +
  # ylab("") +
  guides(fill=guide_legend(title="Type (ICC)", ncol=2),col = FALSE)+
  theme(legend.position="right")+
  theme(legend.text = element_text(size=12))
  


p1_icc



# Load classification by mutation profile (cluster assignment) # CLUSTERS 
cluster_res <- readRDS("../results/euler_memberships_8k_9clusters.rds")

table(cluster_res$clustermembership)
mutation_covariate_data <- read.csv("../data/undivided_binary_only_matrix_6k.csv")

survdata <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")
survdata$clustermembership <- cluster_res$clustermembership
survdata$clustermembership <- as.factor(survdata$clustermembership)
lab <- c("UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2")
levels(survdata$clustermembership) <- lab #LETTERS[1:length(unique(cluster_res$clustermembership))]

survdata <- survdata[!survdata$ELN2022_IPSSM == "IPSSM_NA",]
survdata <- survdata[!survdata$IPSSR_ELN == "IPSSR_NA",]


### SANKEY 
library(networkD3)
# Create a data frame with the two classification vectors
sankeydata <- data.frame(risk = survdata$ELN2022_IPSSM, sample_memberships = survdata$clustermembership)
sankeydata <- na.omit(sankeydata) 
# Create a table with the frequency of each combination
links <- as.data.frame(table(sankeydata))

# Create a nodes data frame
nodes <- data.frame(name = unique(c(as.character(sankeydata$risk), as.character(sankeydata$sample_memberships))))
nodes$id <- 0:(nrow(nodes) - 1)
roworder <-c("ELN2022_adverse","ELN2022_intermediate","ELN2022_favorable","IPSSM_Very-High", "IPSSM_High", "IPSSM_Moderate-High","IPSSM_Moderate-Low","IPSSM_Low", "IPSSM_Very-Low",lab)
nodes <- nodes[match(roworder, nodes$name), ] 


# Map the source and target values in the links data frame to the node IDs
links$Source <- c(match(links$risk, nodes$name) - 1)
links$Target <- c(match(links$sample_memberships, nodes$name) - 1)

links$Target <- match(links$sample_memberships, nodes$name) - 1
links <- links[links$Freq!=0,]

my_color <- 'd3.scaleOrdinal() .domain(["UHR","HR1","NPM1","HR2","INT1","HR3","INT2","LR1","LR2","ELN2022_adverse","ELN2022_intermediate","ELN2022_favorable","IPSSM_Very-High", "IPSSM_High", "IPSSM_Moderate-High","IPSSM_Moderate-Low","IPSSM_Low", "IPSSM_Very-Low"]) .range(["#000000","#772266","#117777","#7c1a29","#114477","#cd9bbc","#88CCAA","#117744","#77AADD","#202020","#774411","#DDAA77","#440154", "#453581","#34618D","#24878E","#88CCAA","#77CC77"])'

#"A", "B", "C","D","E","F","G","H","I"
#"#000000","#774411","#DDAA77","#ed2124","#114477","#CC99BB","#88CCAA","#117744","#77AADD"
#"ELN2022_adverse","ELN2022_intermediate","ELN2022_favorable","IPSSM_Very-High", "IPSSM_High", "IPSSM_Moderate-High","IPSSM_Moderate-Low","IPSSM_Low", "IPSSM_Very-Low"
#"#202020","#774411","#DDAA77","#440154FF", "#453581FF","#34618DFF","#24878EFF","#88CCAA","#77CC77"

library("networkD3")
# Create a Sankey diagram
sankeyPlot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "Source", 
                            Target = "Target", Value = "Freq", NodeID = "name",
                            fontSize = 24, nodeWidth = 25, fontFamily = "sans-serif",
                            colourScale=my_color, iterations = 0, sinksRight = T)

# Display the Sankey diagram
sankeyPlot

# you save it as an html
saveNetwork(sankeyPlot, "../results/sn.html")

library(webshot)
# you convert it as png
webshot("../results/sn.html","./figures/newfig2a.png", vwidth = 900, vheight = 1500)


png(file = "./figures/newfig2b.png", width = 1200, height = 2000, res = 150)
ggarrange(p1_2016, p1_2022, p1_icc, 
          ncol = 1, nrow = 3 , 
          labels = letters[2:4], font.label = c(size=16))
dev.off() 
