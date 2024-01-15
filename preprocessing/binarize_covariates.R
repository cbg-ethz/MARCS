# dichotomize covariates
library("sjmisc")
library("readr")
require("reshape2")
library("plyr")
library(ggplot2)
# load data
input_data <- read.csv("../data/aml_mds_matrix_6k_withCorrectedICC_3.csv",header = T)
colnames(input_data) <- c("ID","ASXL1","ATRX","BCOR","BRAF","CBL","CDKN2A","CREBBP","CUX1",
                          "DNMT3A","ETV6","EZH2","GATA2","GNAS","IDH1","JAK2","KDM6A","KIT",
                          "KRAS","MLL","MPL","MYC","NF1","NPM1","PHF6","PTEN","PTPN11",
                          "RAD21","RUNX1","SF1","SF3B1","SRSF2","STAG2","TET2","TP53","U2AF1",
                          "WT1","ZRSR2","IDH2","NRAS","CEBPA","FLT3","+11","+13","+21",
                          "+22","+8","complex","-12","-17","-18","-20","-3","-5",
                          "-7","-9","inv(16)","inv(3)","-y","t(15;17)","t(6;9)","t(8;21)","t(9;11)",
                          "t(9;22)","t(v;11)","WHO_2016","WHO_2022","ICC","IPSSR_ELN","ELN2022_IPSSM","AGE","SEX","BM_BLASTS","HB","PLT",
                          "WBC","OS","OS_STATUS")
# example output
input_data$AGE[1]
dim(input_data)

# remove columns with NA
# NOTE: there are patients with only few missing values (NAs) but with valid mutations and blood levels.
# Consider if recover them in a second time. 

input_data$ICC[is.na(input_data$ICC)] <- "a"
input_data$WHO_2022[is.na(input_data$WHO_2022)] <- "a"
input_data <- na.omit(input_data)
input_data$ICC[input_data$ICC == "a"] <- NA
input_data$WHO_2022[input_data$WHO_2022 == "a"] <- NA

dim(input_data)


# dichotomize age
binary_covariate_data <- input_data
binary_covariate_data$AGE[grepl("MDS", binary_covariate_data$WHO_2016)] <- as.numeric(cut_number(binary_covariate_data$AGE[grepl("MDS", binary_covariate_data$WHO_2016)], 2))
binary_covariate_data$AGE[grepl("AML", binary_covariate_data$WHO_2016)] <- as.numeric(cut_number(binary_covariate_data$AGE[grepl("AML", binary_covariate_data$WHO_2016)], 2))
binary_covariate_data$AGE[grepl("CMML", binary_covariate_data$WHO_2016)] <- as.numeric(cut_number(binary_covariate_data$AGE[grepl("CMML", binary_covariate_data$WHO_2016)], 2))
#binary_covariate_data$AGE[grepl("aCML", binary_covariate_data$WHO_2016)] <- as.numeric(cut_number(binary_covariate_data$AGE[grepl("aCML", binary_covariate_data$WHO_2016)], 2))
binary_covariate_data$AGE[grepl("aCML", binary_covariate_data$WHO_2016)] <- 1 ## ONE SINGLE CALSE OF aCML here 
binary_covariate_data$AGE[grepl("other", binary_covariate_data$WHO_2016)] <- 1 ## ONE SINGLE CALSE below average here 

# adapt age categories
binary_covariate_data$AGE <- binary_covariate_data$AGE-1

# add sex (1 if male)
binary_covariate_data$SEX <- rep(0,length(input_data$SEX))
binary_covariate_data$SEX[input_data$SEX=="M"] <- 1

# blood
binary_covariate_data$WBC <- as.numeric(cut_number(binary_covariate_data$WBC, 2))-1
binary_covariate_data$PLT <- as.numeric(cut_number(binary_covariate_data$PLT, 2))-1
binary_covariate_data$HB <- as.numeric(cut_number(binary_covariate_data$HB, 2))-1

# bone marrow
binary_covariate_data$BM_BLASTS <- as.numeric(cut(binary_covariate_data$BM_BLASTS, breaks = c(-1,9,19,100)))-1


# put age and sex in last row
col_index <- c(which(colnames(binary_covariate_data)=="AGE"),which(colnames(binary_covariate_data)=="SEX"))
binary_covariate_data <- cbind(binary_covariate_data[, -col_index], binary_covariate_data[, col_index, drop = FALSE])

#write table 
write.csv(binary_covariate_data[,c(2:65,71:74)], file = "../data/binary_only_mutations_blood_6k.csv", quote = F, row.names = binary_covariate_data$ID) ## write matrix
write.csv(binary_covariate_data[,c(66:70,75:78)], file = "../data/binary_only_covariates_OS_6k.csv", quote = F, row.names = binary_covariate_data$ID) ## write matrix
write.csv(binary_covariate_data, file = "../data/undivided_binary_only_matrix_6k.csv", quote = F, row.names=FALSE) # , row.names = binary_covariate_data$ID) 


