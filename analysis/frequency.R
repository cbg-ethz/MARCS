library(dplyr)
df <- survdata[survdata$clustermembership %in% c("HR1","HR2", "HR3"),2:47]

df$WHO_Category <- ifelse(grepl("AML", df$WHO_2016), "AML", 
                          ifelse(grepl("MDS", df$WHO_2016), "MDS", NA))

# Remove rows that do not belong to 'AML' or 'MDS' categories
df <- df %>% filter(!is.na(WHO_Category))

# Initialize a list to store results
results <- list()

# Loop through the gene columns
for (gene in colnames(df)[1:44]) { # Adjust the index if there are more genes
  # Create a binary mutation presence column (0 for no mutation, 1 for mutation)
  df$Mutation_Presence <- ifelse(df[[gene]] > 0, 1, 0)
  
  # Create a 2x2 contingency table
  contingency_table <- table(df$Mutation_Presence, df$WHO_Category)
  
  # Ensure the table is 2x2, if not, add necessary 0s
  if (dim(contingency_table)[1] != 2) {
    contingency_table <- rbind(contingency_table, c(0, 0))
  }
  if (dim(contingency_table)[2] != 2) {
    contingency_table <- cbind(contingency_table, c(0, 0))
  }
  
  # Perform Fisher's exact test
  test <- chisq.test(contingency_table)
  results[[gene]] <- test$p.value
  
  print(paste(colnames(df[gene]),test$p.value))
  print(contingency_table)
  
  # Store the results
  
}

# Create a data frame for the results
results_df <- data.frame(Gene = names(results), P_Value = unlist(results))

# Adjust p-values for multiple testing, if desired, using p.adjust function
results_df$P_Value_Adjusted <- p.adjust(results_df$P_Value, method = "BH")

# Filter for significant results
significant_results <- results_df %>% filter(P_Value < 0.05)

# Print significant genes and their p-values
print(significant_results)
