totdata <- read.csv("/cluster/work/bewi/members/mroncador/diagchallenge/comparison/undivided_binary_only_matrix_8k.csv")
mutation_covariate_data <- read.csv("/cluster/work/bewi/members/mroncador/diagchallenge/comparison/diagCorrected_aml_mds_matrix_8k.csv")

unique(totdata$ID %in% mutation_covariate_data$ID)

totdata <- totdata[order(totdata$ID), ]
mutation_covariate_data <- mutation_covariate_data[order(mutation_covariate_data$ID), ]
totdata$WHO_2016 <- mutation_covariate_data$WHO_2016


### loook for best AIC
library(clustNet)
all_data <- totdata

# remove survival data, classifications and patient id 
index_remove <- c(which(colnames(all_data)=="ID"), which(colnames(all_data)=="OS"),
                  which(colnames(all_data)=="OS_STATUS"), which(colnames(all_data)=="IPSSR_ELN"),
                  which(colnames(all_data)=="WHO_2016"), which(colnames(all_data)=="WHO_2022"),
                  which(colnames(all_data)=="ICC"), which(colnames(all_data)=="ELN2022_IPSSM"))
mut_cov_data <- all_data[,-index_remove]
rownames(mut_cov_data) <- all_data$ID

amlIC <- bestAICsearch(binaryMatrix = mut_cov_data, minK = 5 , maxK = 13)

saveRDS(amlIC, "/cluster/work/bewi/members/mroncador/diagchallenge/comparison/totIC.rds")




