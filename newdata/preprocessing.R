library(readr)
library(dplyr)
library(stringr)
library(devtools)

data <- read.csv(file = "../newdata//MLL_mds_final.csv")
oldata <- read.csv("../data/undivided_binary_only_matrix_6k.csv")

genes <- c(
  "ASXL1", "ATRX", "BCOR", "BRAF", "CBL", "CDKN2A", 
  "CREBBP", "CUX1", "DNMT3A", "ETV6", "EZH2", "GATA2", 
  "GNAS", "IDH1", "JAK2", "KDM6A", "KIT", "KRAS", 
  "MLL", "MPL", "MYC", "NF1", "NPM1", "PHF6", "PTEN", 
  "PTPN11", "RAD21", "RUNX1", "SF1", "SF3B1", "SRSF2", 
  "STAG2", "TET2", "TP53", "U2AF1", "WT1", "ZRSR2", 
  "IDH2", "NRAS", "CEBPA", "FLT3"
)

selectdata <- cbind(data[,1:16], data[,17:56][colnames(data)[17:56] %in% genes], data[,57:61])
selectdata$Cytogenetics.at.sampling <- NULL 
selectdata$Normal <- NULL
selectdata$ANC <- NULL
selectdata$Other <- NULL

colnames(selectdata)[1:13] <- c("ID", "AGE", "IPSSR_ELN", "SEX", "complex", "-5", "-7","-17","-20","+8", "-Y", "OS", "OS_STATUS")
colnames(selectdata)[46:49] <- c("BM_Blasts", "WBC", "HB", "PLT")

## translate months into years
selectdata$OS <- selectdata$OS/12

# clear NAs
selectdata <- selectdata[complete.cases(selectdata),] 

n <- dim(selectdata)[1]

cytogen <- data.frame(
  "t(15;17)" = rep(0, n),
  "t(6;9)" = rep(0, n),
  "t(8;21)" = rep(0, n),
  "t(9;11)" = rep(0, n),
  "t(9;22)" = rep(0, n),
  "t(v;11)" = rep(0, n),
  "inv(16)" = rep(0, n),
  "inv(3)" = rep(0, n)
)

colnames(cytogen) <- c("t(15;17)","t(6;9)","t(8;21)","t(9;11)",
"t(9;22)","t(v;11)","inv(16)","inv(3)")

selectdata <- cbind(selectdata,cytogen)
selectdata$ID <- paste0("CCF_", selectdata$ID)

source("../preprocessing/prepr_functions.R")

selectdata <- classify_aml2016_MLL(selectdata)

selectdata$eln2022_ipssm <- rep(NA, n)

## APPLY IPSSM FOR MDS PATIENTS 
tmpselectdata <- selectdata[selectdata$IPSSR_ELN %in% c("HR-MDS", "LR-MDS"),]
colnames(tmpselectdata)[colnames(tmpselectdata)=="BM_Blasts"] <- "BM_BLAST"
colnames(tmpselectdata)[colnames(tmpselectdata)=="-5"] <- "del5q"
colnames(tmpselectdata)[colnames(tmpselectdata)=="-7"] <- "del7_7q"
colnames(tmpselectdata)[colnames(tmpselectdata)=="-17"] <- "del17_17p"
colnames(tmpselectdata)[colnames(tmpselectdata)=="TP53"] <- "TP53mut"
tmpselectdata$CYTO_IPSSR <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$TP53maxvaf <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$TP53loh <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$MLL_PTD <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$BCORL1 <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$SETBP1 <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$ETNK1 <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$GNB1 <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$ETNK1 <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$PPM1D <- rep(NA, dim(tmpselectdata)[1])
tmpselectdata$PRPF8 <- rep(NA, dim(tmpselectdata)[1])

library(ipssm)
tmpselectdata$BM_BLAST <- as.numeric(tmpselectdata$BM_BLAST)
tmpselectdata$PLT <- as.numeric(tmpselectdata$PLT)
tmpselectdata$HB <- as.numeric(tmpselectdata$HB)
dd.process <- IPSSMprocess(tmpselectdata)
# 3) Calculate IPSS-M
dd.res <- IPSSMmain(dd.process)
# 4) Annotate Results
dd.annot <- IPSSMannotate(dd.res)

myrisk <- dd.annot[,c("ID","IPSSMcat_mean")]
myrisk$IPSSMcat_mean <- paste0("IPSSM_",gsub(" ", "-",x = myrisk$IPSSMcat_mean))

merged_data <- merge(selectdata, myrisk, by = "ID", all.x = TRUE)
selectdata <- merged_data
selectdata$eln2022_ipssm <- selectdata$IPSSMcat_mean
selectdata$IPSSMcat_mean <- NULL

names(selectdata)[c(14:45)] <- tolower(names(selectdata))[c(14:45)]
selectdata <-  calculate_risiko_ELN2022(selectdata)
names(selectdata)[c(14:45)] <- toupper(names(selectdata))[c(14:45)]

names(selectdata)[c(46,59)] <- toupper(names(selectdata))[c(46,59)]
## ADD FEATURES FOR MRC NOT PRESENT, REMOVE AFTER
selectdata$`-12` <- rep(NA, n)
selectdata$`-20` <- rep(NA, n)

selectdata <- classify_icc_who2022(selectdata)
selectdata$`-12` <- NULL 
selectdata$`-20` <- NULL 


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


shared_colnames <- intersect(colnames(input_data),colnames(selectdata))

df1 <- input_data[, shared_colnames]
df2 <- selectdata[, shared_colnames]

totaldf <- rbind(df1,df2)
## Rescue 4 patients with s-AML but blast 0%
totaldf$WHO_2022[totaldf$ID=="CCF_1113"] <- "AML MRC"
totaldf$ICC[totaldf$ID=="CCF_1113"] <- "AML MR-Gen"
totaldf$WHO_2022[totaldf$ID=="CCF_1756"] <- "AML MRC"
totaldf$ICC[totaldf$ID=="CCF_1756"] <- "AML MR-Gen"
totaldf$WHO_2022[totaldf$ID=="CCF_1760"] <- "AML NOS"
totaldf$ICC[totaldf$ID=="CCF_1760"] <- "AML NOS"
totaldf$ICC[totaldf$ID=="CCF_1808"] <- "AML MR-Cyto"

write.csv(totaldf, file = "../data/aml_mds_matrix_8k.csv", quote = F, row.names = F) ## write matrix

# dichotomize covariates
library("sjmisc")
library("readr")
require("reshape2")
library("plyr")
library(ggplot2)

input_data <- na.omit(totaldf)
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

binary_covariate_data$WBC <- as.numeric(binary_covariate_data$WBC)
binary_covariate_data$PLT <- as.numeric(binary_covariate_data$PLT)
binary_covariate_data$HB <- as.numeric(binary_covariate_data$HB)
binary_covariate_data$BM_BLASTS <- as.numeric(binary_covariate_data$BM_BLASTS)
# blood
binary_covariate_data$WBC <- as.numeric(cut_number(binary_covariate_data$WBC, 2))-1
binary_covariate_data$PLT <- as.numeric(cut_number(binary_covariate_data$PLT, 2))-1
binary_covariate_data$HB <- as.numeric(cut_number(binary_covariate_data$HB, 2))-1

# bone marrow
binary_covariate_data$BM_BLASTS <- as.numeric(cut(binary_covariate_data$BM_BLASTS, breaks = c(-1,9,19,100)))-1


# put age and sex in last row
col_index <- c(which(colnames(binary_covariate_data)=="AGE"),which(colnames(binary_covariate_data)=="SEX"))
binary_covariate_data <- cbind(binary_covariate_data[, -col_index], binary_covariate_data[, col_index, drop = FALSE])
binary_covariate_data$IPSSR_ELN <- NULL
binary_covariate_data <- binary_covariate_data[complete_cases(binary_covariate_data),]
write.csv(binary_covariate_data, file = "../data/undivided_binary_only_matrix_8k.csv", quote = F, row.names=FALSE) # , row.names = binary_covariate_data$ID) 

