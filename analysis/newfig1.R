# get bar plots of cluster assignment
rm(list=ls())

library(reshape2)
library(ggplot2)
library(ggpubr)
library(clustNet)

# read data
cluster_results <- readRDS("../results/euler_memberships_8k_9clusters.rds")
mutation_covariate_data <- read.csv("../data/undivided_binary_only_matrix_8k.csv")
data <- read.csv("../data/diagCorrected_aml_mds_matrix_8k_new.csv")

colnames(mutation_covariate_data)[34:46] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
colnames(cluster_results$data)[33:45] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")


for(i in 1:length(cluster_results$DAGs)){
  rownames(cluster_results$DAGs[[i]])[33:45] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
  #rownames(cluster_results$DAGs[[i]])[c(65:70)]<- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")
  
  colnames(cluster_results$DAGs[[i]])[33:45] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
  #colnames(cluster_results$DAGs[[i]])[c(65:70)]<- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")
}

# plot graphs with constant node size
#graphClust::plot_clusters(cluster_results)

# and with entropy node size
all_data <- mutation_covariate_data
index_remove <- c(which(colnames(all_data)=="ID"), which(colnames(all_data)=="OS"), which(colnames(all_data)=="OS_STATUS"), which(colnames(all_data)=="IPSSR_ELN"), which(colnames(all_data)=="WHO_2016"),which(colnames(all_data)=="WHO_2022"),which(colnames(all_data)=="ICC"),which(colnames(all_data)=="ELN2022_IPSSM"))
mut_cov_data <- all_data[,-index_remove]
node_colours <- c(rep("#92a8d1",32),rep("#064273",13) , rep("#ff727c",4), rep("#708090",2))


netplot <- clustNet::plot_clusters(cluster_results,  directed = F, node_colours = node_colours)

plot_clusters2 <- function (cluster_results, node_colours = "#fdae61", scale_entropy = FALSE, 
          directed = TRUE) 
{
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package \"ggpubr\" must be installed to use this function.", 
         call. = FALSE)
  }
  data <- as.data.frame(cluster_results$data)
  node_size <- NULL
  if (!is.null(data)) {
    k_clust <- length(cluster_results$DAGs)
    total_entropies <- sapply(data, function(x) entropy(x))
    cluster_dim <- c()
    diff_entropies <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    binary_frequency <- matrix(NA, nrow = k_clust, ncol = ncol(data))
    for (nn in 1:k_clust) {
      cluster_data <- data[cluster_results$clustermembership == 
                             nn, ]
      cluster_dim[nn] <- dim(cluster_data)[1]
      cluster_entropies <- sapply(cluster_data, function(x) entropy(x))
      diff_entropies[nn, ] <- cluster_entropies - total_entropies
      binary_frequency[nn, ] <- sapply(cluster_data, function(x) sum(x))
    }
    binary_frequency <- sweep(binary_frequency, 1, cluster_dim, 
                              FUN = "/")
    binary_frequency_temp <- binary_frequency
    max_col_values <- apply(binary_frequency_temp, 2, max)
    node_size_percentage_frequ_temp <- sweep(binary_frequency_temp, 
                                             2, max_col_values, FUN = "/")
    node_size <- 1.5 + node_size_percentage_frequ_temp * 
      4
    if (scale_entropy) {
      node_size_percentage <- 1 - (diff_entropies[, ] + 
                                     abs(min(diff_entropies[, ])))/(max(diff_entropies[, 
                                     ]) + abs(min(diff_entropies[, ])))
      node_size <- 1.5 + node_size_percentage * 6
    }
  }
  p_list <- list()
  k_clust <- length(cluster_results$DAGs)
  for (ii in 1:k_clust) {
    my_DAG <- cluster_results$DAGs[ii][[1]]
    p_list[[ii]] <- nice_DAG_plot(my_DAG, print_direct = FALSE, 
                                  node_size = node_size[ii, ], node_colours = node_colours, 
                                  directed = directed)
  }
  lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")
  ggpubr::ggarrange(plotlist = p_list, labels = paste(lab)) #"Cluster", LETTERS[1:k_clust],"-",
}

library(entropy)
netplot2 <- plot_clusters2(cluster_results,  directed = F, node_colours = node_colours)

netplot2

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

#colnames(mutation_covariate_data)[43:65] <- c("+11","+13","+21", "+22", "+8", "complex", "-12", "del(17p)/-17","-18", "-20", "-3", "-5/-5q", "-7/-7q", "-9", "inv(16)","inv(3)",
#                                              "-Y", "t(15;17)", "t(6;9)", "t(8;21)",  "t(9;11)",  "t(9;22)", "t(v;11)")

colnames(mutation_covariate_data)[34:46] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")

#colnames(mutation_covariate_data)[c(71:74,77,78)] <- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")

#colnames(cluster_results$data)[42:64] <- c("+11","+13","+21", "+22", "+8", "complex", "-12", "del(17p)/-17","-18", "-20", "-3", "-5/-5q", "-7/-7q", "-9", "inv(16)","inv(3)",
# "-Y", "t(15;17)", "t(6;9)", "t(8;21)",  "t(9;11)",  "t(9;22)", "t(v;11)")
#colnames(cluster_results$data)[c(65:70)] <- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")

colnames(cluster_results$data)[33:45] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")


for(i in 1:length(cluster_results$DAGs)){
  rownames(cluster_results$DAGs[[i]])[33:45] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
  #rownames(cluster_results$DAGs[[i]])[c(65:70)]<- c("BM Blasts", "Hb", "Plt", "Wbc","Age","Sex")
  
  colnames(cluster_results$DAGs[[i]])[33:45] <- c("+8", "complex", "del(17p)/-17", "del(5q)/-5", "del(7q)/-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
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
# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Cluster"
levels(melted_data_sum$Cluster) <- LETTERS[1:length(unique(cluster_results$clustermembership))]

lab <- c("UHR","HR1","NPM1","HR2","INT2","HR3","INT1","LR1","LR2")

melted_data_sum$Cluster[melted_data_sum$Cluster=="1"] <- lab[1]
melted_data_sum$Cluster[melted_data_sum$Cluster=="2"] <- lab[2]
melted_data_sum$Cluster[melted_data_sum$Cluster=="3"] <- lab[3]
melted_data_sum$Cluster[melted_data_sum$Cluster=="4"] <- lab[4]
melted_data_sum$Cluster[melted_data_sum$Cluster=="5"] <- lab[5]
melted_data_sum$Cluster[melted_data_sum$Cluster=="6"] <- lab[6]
melted_data_sum$Cluster[melted_data_sum$Cluster=="7"] <- lab[7]
melted_data_sum$Cluster[melted_data_sum$Cluster=="8"] <- lab[8]
melted_data_sum$Cluster[melted_data_sum$Cluster=="9"] <- lab[9]


#colours_clusters <- c("#000000","#772266","#117777","#7c1a29","#114477","#cd9bbc","#88CCAA","#117744","#77AADD")
colours_clusters <- c("#2f0000","#9B2226","#94D2BD","#BB3E03","#E9D8A6","#CA6702","#EE9B00","#0A9396","#005F73")




# Stacked barplot (sum of mutations)
barplot_sum0 <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("UHR"="#1f0000","HR1"="#9B2226","NPM1"="#94D2BD","HR2"="#BB3E03","INT2"="#E9D8A6","HR3"="#CA6702","INT1"= "#EE9B00","LR1"="#0A9396","LR2"="#005F73"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10))+
  theme(axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 12))+
  ggtitle(label = paste("Genetics")) + theme(plot.title = element_text(color="black", size=13, face="bold"))+
  theme(legend.text=element_text(size=10), legend.title = element_text(size = 10))


  
pdf(file = "./figures/newfig1.pdf", width = 14, height = 17)
ggarrange(netplot2, barplot_sum0, 
          ncol = 1, nrow = 2 , 
          labels = letters[1:6], 
          heights = c(1, .32))
dev.off() 

