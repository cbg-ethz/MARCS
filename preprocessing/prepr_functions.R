
parse_cytogenetics <- function(df) {
  if (!"CYTOGENETICS" %in% colnames(df)) {
    stop("The input data frame must have a 'CYTOGENETICS' column.")
  }
  
  # define alterations of interest
  alterations <- c("t\\(15;17\\)\\(q24;q21\\)", "t\\(8;21\\)\\(q22;q22\\)", "inv\\(16\\)\\(p13.1q22\\)",
                   "t\\(9;11\\)\\(p21;q23\\)", "t\\(6;9\\)\\(p22;q34\\)", "t\\(9;22\\)\\(q34;q11.2\\)", "\\+11", "\\+13",
                   "\\+21", "\\+22", "\\+8", "-12", "-17", "-18", "-20", "-3", "-4", "-9", "-Y",
                   "t\\(.*;11\\)", "t\\(12;21\\)\\(q24;q22\\)", "t\\(3;21\\)\\(q26;q22\\)",
                   "del\\(17\\)\\(q21\\)", "t\\(6;11\\)\\(q27;q23\\)")
  
  # alterations in a human readable format
  readable_alterations <- c("t(15;17)", "t(8;21)", "inv(16)", "t(9;11)",
                            "t(6;9)", "t(9;22)",
                            "+11", "+13", "+21", "+22", "+8", "-12", "-17", "-18", "-20", "-3", 
                            "-4", "-9", "-Y", "t(v;11)", "t(12;21)",  "t(3;21)",
                            "del(17q)", "t(6;11)")
  
  # sanitize 'CYTOGENETICS' column: remove whitespace
  df$CYTOGENETICS <- gsub("\\s", "", df$CYTOGENETICS)
  
  # create new columns for each alteration
  for (i in seq_along(alterations)) {
    colname <- readable_alterations[i] #gsub("\\(|\\)|;", "_", readable_alterations[i])
    df[[colname]] <- as.integer(grepl(alterations[i], df$CYTOGENETICS))
  }
  
  # merge information for del(5q) and -5, and del(7q) and -7
  df$`-5` <- as.integer(grepl("del\\(5q\\)", df$CYTOGENETICS) | grepl("del\\(5\\)\\(q", df$CYTOGENETICS) | grepl("-5", df$CYTOGENETICS))
  df$`-7` <- as.integer(grepl("del\\(7q\\)", df$CYTOGENETICS) | grepl("del\\(7\\)\\(q", df$CYTOGENETICS) | grepl("-7", df$CYTOGENETICS))
  df$`del(12p)` <- as.integer(grepl("del\\(12p\\)", df$CYTOGENETICS) | grepl("del\\(12\\)\\(p", df$CYTOGENETICS) | grepl("-12", df$CYTOGENETICS))
  df$`inv(3)` <- as.integer(grepl("inv\\(3\\)\\(q21q26\\)", df$CYTOGENETICS) | grepl("t\\(3;3\\)\\(q21;q26\\)", df$CYTOGENETICS) )
  
  
  
  # create 'complex' column
  df$complex <- apply(df, 1, function(x) {
    cytogenetics <- strsplit(x["CYTOGENETICS"], split = "[/,]")[[1]]
    return(as.integer(length(cytogenetics) > 5))
  })
  
  # remove 'CYTOGENETICS' column
  df$CYTOGENETICS <- NULL
  
  return(df)
}




#### FUNCTION FOR CLASSIFICATION WHO 2016
library(dplyr)

classify_aml2016_nejm <- function(df) {
  df <- df %>%
    mutate(WHO_2016 = case_when(
      `t(8;21)` == 1 ~ "AML t(8;21)",
      `inv(16)` == 1 ~ "AML inv(16)",
      `t(15;17)` == 1 ~ "AML t(15;17)",
      BM_Blasts > 20 & `t(6;9)` == 1 ~ "AML t(6;9)",
      BM_Blasts > 20 & `inv(3)` == 1 ~ "AML inv(3)",
      BM_Blasts > 20 & `t(9;11)` == 1 ~ "AML t(9;11)",
      BM_Blasts > 20 & RUNX1 == 1 ~ "AML RUNX1",
      BM_Blasts > 20 & NPM1 == 1 ~ "AML NPM1",
      BM_Blasts > 20 & CEBPA == 2 ~ "AML CEBPAbi",
      BM_Blasts > 20 & (`-5` == 1 | `-7` == 1 | complex == 1 ) ~ "AML MRC",
      TRUE ~ "AML NOS"
    ))
  
  return(df)
}

classify_aml2016_MLL <- function(df) {
  df <- df %>%
    mutate(WHO_2016 = case_when(
      `t(8;21)` == 1 ~ "AML t(8;21)",
      `inv(16)` == 1 ~ "AML inv(16)",
      `t(15;17)` == 1 ~ "AML t(15;17)",
      BM_Blasts >= 20 & `t(6;9)` == 1 ~ "AML t(6;9)",
      BM_Blasts >= 20 & `inv(3)` == 1 ~ "AML inv(3)",
      BM_Blasts >= 20 & `t(9;11)` == 1 ~ "AML t(9;11)",
      BM_Blasts >= 20 & RUNX1 == 1 ~ "AML RUNX1",
      BM_Blasts >= 20 & NPM1 == 1 ~ "AML NPM1",
      BM_Blasts >= 20 & CEBPA == 2 ~ "AML CEBPAbi",
      BM_Blasts >= 20 & (`-5` == 1 | `-7` == 1 | complex == 1 ) ~ "AML MRC",
      BM_Blasts >= 20 ~ "AML NOS",
      BM_Blasts < 20 ~ IPSSR_ELN
    ))
  
  return(df)
}




classify_aml2016_ncri <- function(df) {
  df <- df %>%
    mutate(WHO_2016 = case_when(
      `t(8;21)` == 1 ~ "AML t(8;21)",
      `inv(16)` == 1 ~ "AML inv(16)",
      `t(15;17)` == 1 ~ "AML t(15;17)",
      bm_blasts > 20 & `t(6;9)` == 1 ~ "AML t(6;9)",
      bm_blasts > 20 & `inv(3)` == 1 ~ "AML inv(3)",
      bm_blasts > 20 & `t(9;11)` == 1 ~ "AML t(9;11)",
      bm_blasts > 20 & RUNX1 == 1 ~ "AML RUNX1",
      bm_blasts > 20 & NPM1 == 1 ~ "AML NPM1",
      bm_blasts > 20 & CEBPA == 2 ~ "AML CEBPAbi",
      bm_blasts > 20 & (`-5` == 1 | `-7` == 1 | complex == 1 | `t(3;21)` == 1 | `t(1;3)` == 1 | `t(5;12)` == 1 | `t(5;17)` == 1 | `t(3;5)` == 1) ~ "AML MRC",
      TRUE ~ "AML NOS"
    ))
  
  return(df)
}



classify_icc_who2022 <- function(df) {
  
  df <- df %>%
    mutate(
      WHO_2022 = case_when(
        str_detect(WHO_2016, "MDS/MPN-U") ~ "MDS/MPN-U",
        str_detect(WHO_2016, "MDS/MPN-RS-T") ~ "MDS/MPN-RS-T",
        str_detect(WHO_2016, "aCML") ~ "aCML",
        `t(8;21)` == 1 ~ "AML t(8;21)",
        `inv(16)` == 1 ~ "AML inv(16)",
        `t(15;17)` == 1 ~ "AML t(15;17)",
        `t(6;9)` == 1 ~ "AML t(6;9)",
        `inv(3)` == 1 ~ "AML inv(3)",
        `t(9;11)` == 1 ~ "AML t(v;11)",
        `t(v;11)` == 1 ~ "AML t(v;11)",
        NPM1 == 1 ~ "AML NPM1",
        BM_BLASTS > 20 & CEBPA == 2 ~ "AML CEBPAbi",
        BM_BLASTS > 20 & (`-5` == 1 | `-7` == 1 | `-17` == 1 | `+8` == 1 | `-12` == 1 | complex == 1 |
                            ASXL1 == 1 | BCOR == 1 | EZH2 == 1 | SF3B1 == 1 | SRSF2 == 1 | STAG2 == 1 | U2AF1 == 1 | ZRSR2 == 1) ~ "AML MRC",
        BM_BLASTS >= 20 ~ "AML NOS",
        str_detect(WHO_2016, "CMML") ~ "CMML",
        BM_BLASTS < 5 & `-5` == 1 ~ "MDS-5q",
        BM_BLASTS < 20 & SF3B1 == 1 & `-5` == 0 & `-7` == 0 & complex == 0 ~ "MDS-SF3B1",
        BM_BLASTS < 20 & TP53 == 2 ~ "MDS-biTP53",
        str_detect(WHO_2016, "MDS") & BM_BLASTS < 5 ~ "MDS-LB",
        BM_BLASTS >= 5 & BM_BLASTS < 10 ~ "MDS-IB1",
        BM_BLASTS >= 10 & BM_BLASTS < 20 ~ "MDS-IB2"
      ),
      ICC = case_when(
        str_detect(WHO_2016, "MDS/MPN-U") ~ "MDS/MPN-U",
        str_detect(WHO_2016, "MDS/MPN-RS-T") ~ "MDS/MPN-RS-T",
        str_detect(WHO_2016, "aCML") ~ "aCML",
        `t(8;21)` == 1 & BM_BLASTS > 10 ~ "AML t(8;21)",
        `inv(16)` == 1 & BM_BLASTS > 10 ~ "AML inv(16)",
        `t(15;17)` == 1 & BM_BLASTS > 10 ~ "AML t(15;17)",
        `t(6;9)` == 1 & BM_BLASTS > 10 ~ "AML t(6;9)",
        `inv(3)` == 1 & BM_BLASTS > 10 ~ "AML inv(3)",
        `t(9;11)` == 1 ~ "AML t(9;11)",
        `t(9;22)` == 1 & BM_BLASTS > 20 ~ "AML t(9;22)",
        NPM1 == 1 & BM_BLASTS > 10 ~ "AML NPM1",
        CEBPA == 2 & BM_BLASTS > 10 ~ "AML CEBPAbi",
        TP53 == 1 & BM_BLASTS > 20 ~ "AML TP53",
        BM_BLASTS >= 20 & (ASXL1 == 1 | BCOR == 1 | EZH2 == 1 | RUNX1 == 1| SF3B1 == 1 | SRSF2 == 1 | STAG2 == 1 | U2AF1 == 1 | ZRSR2 == 1) ~ "AML MR-Gen",
        BM_BLASTS >= 20 & (`-5` == 1 | `-7` == 1 | `-17` == 1 | `+8` == 1 | `-12` == 1 | `-20` == 1 | complex == 1) ~ "AML MR-Cyto",
        BM_BLASTS >= 20 ~ "AML NOS",
        TP53 == 1 & BM_BLASTS >= 10 & BM_BLASTS < 20 ~ "MDS/AML TP53",
        BM_BLASTS >= 10 & BM_BLASTS < 20 & (`-5` == 1 | `-7` == 1 | `-17` == 1 | `+8` == 1 | `-12` == 1 | `-20` == 1 | complex == 1) ~ "MDS/AML MR-Cyto",
        BM_BLASTS >= 10 & BM_BLASTS < 20 & (ASXL1 == 1 | BCOR == 1 | EZH2 == 1 | RUNX1 == 1| SF3B1 == 1 | SRSF2 == 1 | STAG2 == 1 | U2AF1 == 1 | ZRSR2 == 1) ~ "MDS/AML MR-Gen",
       `-5` == 1 & BM_BLASTS < 10 ~ "MDS-5q",
        SF3B1 == 1 & TP53 < 2 & RUNX1 == 0 & `-5` == 0 & `-7` == 0 & complex == 0 & BM_BLASTS < 20 ~ "MDS-SF3B1",
        (TP53 == 2 | TP53 == 1 & complex==1) & BM_BLASTS < 10 ~ "MDS-TP53",
        str_detect(WHO_2016, "MDS") & BM_BLASTS < 5 ~ "MDS-LB",
        BM_BLASTS >= 5 & BM_BLASTS < 10 ~ "MDS-EB",
        NPM1 == 0 & BM_BLASTS >= 10 & BM_BLASTS < 20 ~ "MDS/AML NOS"
      )
    )
  
  return(df)
}



calculate_risiko_ELN2022 <- function(dataframe) {
  
  favorable_genes <- c("t(8;21)", "inv(16)", "npm1")
  intermediate_genes <- c("flt3", "t(9;11)")
  adverse_genes <- c("t(6;9)", "t(v;11)", "t(9;22)", "inv(3)", "-5",
                     "i(17)", "-7", "complex", "asxl1", "bcor",
                     "ezh2", "runx1", "sf3b1", "srsf2", "stag2", "u2af1",
                     "zrsr2", "tp53")
  
  for (i in 1:nrow(dataframe)) {
    # Skip if there is already a value in eln2022_ipssm for this row
    if (!is.na(dataframe$eln2022_ipssm[i])) next
    
    gene_mutations <- colnames(dataframe)[which(dataframe[i,] == 1)]
    
    if (any(gene_mutations %in% adverse_genes)) {
    dataframe$eln2022_ipssm[i] <- "ELN2022_adverse"
    } else if (dataframe$flt3[i]==1 & dataframe$npm1[i]==1) {
      dataframe$eln2022_ipssm[i] <- "ELN2022_intermediate"
    } else if (dataframe$flt3[i]==1 & (!any(gene_mutations %in% adverse_genes))) {
      dataframe$eln2022_ipssm[i] <- "ELN2022_intermediate"
    } else if (dataframe$`t(9;11)`[i]==1) {
      dataframe$eln2022_ipssm[i] <- "ELN2022_intermediate"
    } else if (dataframe$`t(8;21)`[i]==1) {
      dataframe$eln2022_ipssm[i] <- "ELN2022_favorable"
    } else if (dataframe$`inv(16)`[i]==2) {
      dataframe$eln2022_ipssm[i] <- "ELN2022_favorable"
    } else if (dataframe$npm1[i]==1 & dataframe$flt3[i]==0) {
      dataframe$eln2022_ipssm[i] <- "ELN2022_favorable"
    } else if (dataframe$cebpa[i]>=1) {
      dataframe$eln2022_ipssm[i] <- "ELN2022_favorable"
    } else {
      dataframe$eln2022_ipssm[i] <- "ELN2022_intermediate"
    }
  }
  
  return(dataframe)
}




