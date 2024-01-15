library(ggplot2)
library(dplyr)
library(tidyr)


data <- read.csv(file = "../data/aml_mds_matrix_8k.csv")
data2 <-  read.csv("../data/diagCorrected_aml_mds_matrix_8k.csv")
data <- data[data$ID %in% data2$ID,]
colnames(data)[34:46] <- c("+8", "complex", "-17", "-5", "-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")
data$CEBPA[data$CEBPA==2] <- 1

library(dplyr)

# Assuming your data frame is named data
data <- data %>%
  mutate(Disease = ifelse(grepl("AML", WHO_2016), "AML", 
                          ifelse(grepl("MDS", WHO_2016), "MDS", NA)))

subdata <- data[,c(2:46,60)]


mutation_sum <- aggregate(subdata[,1:45], by = list(subdata$Disease), FUN = sum)

mutation_sum <- as.data.frame(mutation_sum)
mutation_sum$`inv(16)` <- NULL 
mutation_sum$`inv(3)` <- NULL 
mutation_sum$`t(15;17)` <- NULL 
mutation_sum$`t(6;9)` <- NULL 
mutation_sum$`t(8;21)` <- NULL 
mutation_sum$`t(9;11)` <- NULL 
mutation_sum$`t(9;22)` <- NULL 
mutation_sum$`t(v;11)` <- NULL 

# Reshape the data for plotting
melted_data <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data)[1] <- "Cluster"

# Reshape the data for stacked barplot (sum of mutations)
melted_data_sum <- melt(mutation_sum, id.vars = "Group.1")
colnames(melted_data_sum)[1] <- "Disease"

for(i in seq(1,dim(melted_data_sum)[1],2)){
  melted_data_sum$tot[i] <- as.numeric(melted_data_sum[i,3]+melted_data_sum[i+1,3])
  melted_data_sum$tot[i+1] <- as.numeric(melted_data_sum[i,3]+melted_data_sum[i+1,3]-1)
}

melted_data_sum <- melted_data_sum[order(-melted_data_sum$tot),]

melted_data_sum$variable <- factor(melted_data_sum$variable, unique(melted_data_sum$variable))

# Stacked barplot (sum of mutations)
barplot_sum <- ggplot(melted_data_sum, aes(x = variable, y = value, fill = Disease)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("AML"="#440154FF","MDS"="#35B779FF"))+
  xlab("Genes") +
  ylab("Sum of mutations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =10))+
  theme(axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 10))+
  ggtitle(label = paste("Mutations in MDS/AML cohort (n=7480)")) + theme(plot.title = element_text(color="black", size=12, face="bold"))+
  theme(legend.text=element_text(size=10), legend.title = element_text(size = 10))+
  guides(fill = guide_legend(title = "Disease (WHO2016)"))
barplot_sum

png(filename = "../results/genetic_landscape.png", width = 1200, height = 500, res = 120)
barplot_sum
dev.off()


subdata2 <- data[,c(54:57,60)]
subdata2$Disease[is.na(subdata2$Disease)] <- "other"
## REMOVE ERRORS IN TYPING 
subdata2$HB[subdata2$HB>20] <- subdata2$HB[subdata2$HB>20]/10
  
mean(subdata2$BM_BLASTS[subdata2$Disease=="AML"])
mean(subdata2$BM_BLASTS[subdata2$Disease=="MDS"])
mean(subdata2$HB[subdata2$Disease=="AML"])
mean(subdata2$HB[subdata2$Disease=="MDS"])
mean(subdata2$WBC[subdata2$Disease=="AML"])
mean(subdata2$WBC[subdata2$Disease=="MDS"])
mean(subdata2$PLT[subdata2$Disease=="AML"])
mean(subdata2$PLT[subdata2$Disease=="MDS"])

t.test(subdata2$BM_BLASTS[subdata2$Disease=="AML"], subdata2$BM_BLASTS[subdata2$Disease=="MDS"])
t.test(subdata2$HB[subdata2$Disease=="AML"], subdata2$HB[subdata2$Disease=="MDS"])
t.test(subdata2$WBC[subdata2$Disease=="AML"], subdata2$WBC[subdata2$Disease=="MDS"])
t.test(subdata2$PLT[subdata2$Disease=="AML"], subdata2$PLT[subdata2$Disease=="MDS"])


# Function to create individual histograms with thin lines
plot_individual_column <- function(col_name) {
  ggplot(subdata2, aes_string(x = col_name, color = "Disease")) +
    geom_freqpoly(binwidth = 1, size = 1) +
    labs(x = col_name, y = "Frequency", color = "Disease") +
    scale_color_discrete(name = "Disease") +
    theme_minimal()
}

# List of column names for which you want to create histograms
columns_to_plot <- c("BM_BLASTS", "HB", "WBC", "PLT")

# Create and arrange individual plots in a grid
plots_list <- lapply(columns_to_plot, plot_individual_column)

png(filename = "../results/bloodvalues_landscape.png", width = 1200, height = 500, res = 120)
gridExtra::grid.arrange(grobs = plots_list, ncol = 2)
dev.off()





