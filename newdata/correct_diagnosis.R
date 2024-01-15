data <- read.csv(file = "../data/aml_mds_matrix_8k.csv")
data2 <-  read.csv("../data/undivided_binary_only_matrix_8k.csv")
data <- data[data$ID %in% data2$ID,]
colnames(data)[34:46] <- c("+8", "complex", "-17", "-5", "-7","inv(16)","inv(3)","t(15;17)","t(6;9)", "t(8;21)","t(9;11)","t(9;22)","t(v;11)")


asd <- data[data$IPSSR_ELN %in% c("s-AML","LR-MDS","HR-MDS"),] %>%
  mutate(WHO_2016 = case_when(
    `t(8;21)` == 1 ~ "AML t(8;21)",
    `inv(16)` == 1 ~ "AML inv(16)",
    `t(15;17)` == 1 ~ "AML t(15;17)",
    BM_BLASTS > 20 & `t(6;9)` == 1 ~ "AML t(6;9)",
    BM_BLASTS > 20 & `inv(3)` == 1 ~ "AML inv(3)",
    BM_BLASTS > 20 & `t(9;11)` == 1 ~ "AML t(9;11)",
    BM_BLASTS > 20 & RUNX1 == 1 ~ "AML RUNX1",
    BM_BLASTS > 20 & NPM1 == 1 ~ "AML NPM1",
    BM_BLASTS > 20 & CEBPA == 2 ~ "AML CEBPAbi",
    BM_BLASTS > 20 & (`-5` == 1 | `-7` == 1 | complex == 1 ) ~ "AML MRC",
    BM_BLASTS > 10 ~ "MDS-EB2",
    BM_BLASTS > 5 ~ "MDS-EB1",
    IPSSR_ELN == "LR-MDS" & (`-5` == 1) ~ "MDS-del5q",
    IPSSR_ELN=="s-AML" ~ "AML MRC",
    TRUE ~ IPSSR_ELN
  ))

ids <- data$ID %in% asd$ID
data$WHO_2016[ids] <- asd$WHO_2016

signif(table(data$WHO_2016)*100/7480,digits = 4)

### NOW CORRECT FOR THE CORRECTED ICC 
#iccinfo <- read.csv(file = "../data/aml_mds_matrix_6k_withCorrectedICC_2.csv")
iccinfo <- read.csv(file = "../data/aml_mds_matrix_6k_withCorrectedICC_3.csv")
iccinfo <- iccinfo[iccinfo$ID %in% data2$ID,]
sharedid <- data$ID[data$ID %in% iccinfo$ID]
data$ICC[data$ID %in% sharedid] <- iccinfo$ICC
data$ICC[data$WHO_2016=="CMML"] <- "CMML"
data$WHO_2022[data$ID %in% sharedid] <- iccinfo$WHO_2022

write.csv(data, file = "../data/diagCorrected_aml_mds_matrix_8k_new.csv", quote = F, row.names=FALSE) # , row.names = binary_covariate_data$ID) 
