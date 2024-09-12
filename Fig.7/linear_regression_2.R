# Clear environment
rm(list = ls())
options(stringsAsFactors=FALSE)

# Set path and working directory
path <- "/share1/kitani/data_from_first/elife_2017"
setwd(path)

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(tibble)
library(lmtest)

# Define directories
out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"

# Load datasets
data <- read.csv("gene_expression_matrix_2016-03-03/columns-samples.csv")
exp <- read.csv("gene_expression_matrix_2016-03-03/fpkm_table_unnormalized.csv", row.names = 1)
gene <- read.csv("gene_expression_matrix_2016-03-03/rows-genes.csv")
protein_df <- read.csv("ProteinAndPathologyQuantifications.csv")
meta_df <- read.csv("DonorInformation.csv")

# Clean column names
names(exp) <- gsub("X", "", names(exp))

# Extract gene expression for GPNMB
gpnmb_id <- gene[gene$gene_symbol == "GPNMB", ]$gene_id
exp_GPNMB <- exp[as.character(gpnmb_id), ]
exp_GPNMB <- exp_GPNMB[ , -which(names(exp_GPNMB) == "496100436")]  # Remove specific column
exp_GPNMB_t <- as.data.frame(t(exp_GPNMB))  # Transpose
exp_GPNMB_t$rnaseq_profile_id <- row.names(exp_GPNMB_t)
exp_GPNMB_t <- merge(exp_GPNMB_t, data, by = "rnaseq_profile_id")
names(exp_GPNMB_t)[2] <- "GPNMB"

# Merge datasets
combined <- merge(exp_GPNMB_t, meta_df, by = "donor_id")
combined_2 <- merge(combined, protein_df, by = c("donor_id", "structure_id"))

# Build linear regression model
model <- lm(ihc_tau2_ffpe ~ ihc_a_beta_ffpe * GPNMB, data = combined_2)
summary(model)

# Loop through structures and diagnoses for model fitting
unique_structures <- unique(combined_2$structure_acronym.x)
unique_diagnoses <- unique(combined_2$act_demented)

for (struct in unique_structures) {
  for (diag in unique_diagnoses) {
    
    subset_data <- subset(combined_2, structure_acronym.x == struct & act_demented == diag)
    
    # Ensure sufficient data
    if (nrow(subset_data) > 1) {
      model <- lm(ihc_tau2_ffpe ~ ihc_a_beta_ffpe * GPNMB, data = subset_data)
      cat("Summary for structure:", struct, "and diagnosis:", diag, "\n")
      print(summary(model))
      cat("\n\n")
    } else {
      cat("Insufficient data for structure:", struct, "and diagnosis:", diag, "\n\n")
    }
  }
}

# Function to extract summary statistics from models
extract_summary_stats <- function(df, struct, diag) {
  
  # Filter by structure and diagnosis
  subset_data <- subset(df, structure_acronym.x == struct & act_demented == diag)
  
  # Ensure sufficient data
  if (nrow(subset_data) <= 1) {
    return(data.frame())
  }
  
  # Build the model
  model <- lm(ihc_tau2_ffpe ~ ihc_a_beta_ffpe * GPNMB, data = subset_data)
  summary_result <- summary(model)
  coefs <- summary_result$coefficients
  
  # Extract relevant statistics
  data.frame(
    structure_acronym = struct,
    act_demented = diag,
    intercept_estimate = coefs["(Intercept)", "Estimate"],
    intercept_stderror = coefs["(Intercept)", "Std. Error"],
    a_beta_estimate = coefs["ihc_a_beta_ffpe", "Estimate"],
    a_beta_stderror = coefs["ihc_a_beta_ffpe", "Std. Error"],
    gpnmb_estimate = coefs["GPNMB", "Estimate"],
    gpnmb_stderror = coefs["GPNMB", "Std. Error"],
    interaction_estimate = ifelse("ihc_a_beta_ffpe:GPNMB" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:GPNMB", "Estimate"], NA),
    interaction_stderror = ifelse("ihc_a_beta_ffpe:GPNMB" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:GPNMB", "Std. Error"], NA),
    r_squared = summary_result$r.squared,
    adj_r_squared = summary_result$adj.r.squared,
    p_value_intercept = coefs["(Intercept)", "Pr(>|t|)"],
    p_value_a_beta = coefs["ihc_a_beta_ffpe", "Pr(>|t|)"],
    p_value_gpnmb = coefs["GPNMB", "Pr(>|t|)"],
    p_value_interaction = ifelse("ihc_a_beta_ffpe:GPNMB" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:GPNMB", "Pr(>|t|)"], NA)
  )
}

# Loop through structures and diagnoses to extract stats
results_list <- list()
for (struct in unique_structures) {
  for (diag in unique_diagnoses) {
    result <- extract_summary_stats(combined_2, struct, diag)
    if (nrow(result) > 0) {
      results_list[[length(results_list) + 1]] <- result
    } else {
      cat("Insufficient data for structure:", struct, "and diagnosis:", diag, "\n\n")
    }
  }
}

# Combine all results into a single data frame and write to CSV
all_stats <- do.call(rbind, results_list)
write.csv(all_stats, "out.dir.kitani/plot/lm_all_stats_GPNMB.csv", row.names = FALSE, quote = F)


#################### AT8 Analysis ###################

# Function to extract summary statistics for AT8
extract_summary_stats_at8 <- function(df, struct, diag) {
  
  # Filter by structure and diagnosis
  subset_data <- subset(df, structure_acronym.x == struct & act_demented == diag)
  
  # Ensure sufficient data
  if (nrow(subset_data) <= 1) {
    return(data.frame())
  }
  
  # Build the model
  model <- lm(ihc_at8_ffpe ~ ihc_a_beta_ffpe * GPNMB, data = subset_data)
  summary_result <- summary(model)
  coefs <- summary_result$coefficients
  
  # Extract relevant statistics
  data.frame(
    structure_acronym = struct,
    act_demented = diag,
    intercept_estimate = coefs["(Intercept)", "Estimate"],
    intercept_stderror = coefs["(Intercept)", "Std. Error"],
    a_beta_estimate = coefs["ihc_a_beta_ffpe", "Estimate"],
    a_beta_stderror = coefs["ihc_a_beta_ffpe", "Std. Error"],
    gpnmb_estimate = coefs["GPNMB", "Estimate"],
    gpnmb_stderror = coefs["GPNMB", "Std. Error"],
    interaction_estimate = ifelse("ihc_a_beta_ffpe:GPNMB" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:GPNMB", "Estimate"], NA),
    interaction_stderror = ifelse("ihc_a_beta_ffpe:GPNMB" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:GPNMB", "Std. Error"], NA),
    r_squared = summary_result$r.squared,
    adj_r_squared = summary_result$adj.r.squared,
    p_value_intercept = coefs["(Intercept)", "Pr(>|t|)"],
    p_value_a_beta = coefs["ihc_a_beta_ffpe", "Pr(>|t|)"],
    p_value_gpnmb = coefs["GPNMB", "Pr(>|t|)"],
    p_value_interaction = ifelse("ihc_a_beta_ffpe:GPNMB" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:GPNMB", "Pr(>|t|)"], NA)
  )
}

# Loop through structures and diagnoses to extract stats for AT8
results_list_at8 <- list()
for (struct in unique_structures) {
  for (diag in unique_diagnoses) {
    result_at8 <- extract_summary_stats_at8(combined_2, struct, diag)
    if (nrow(result_at8) > 0) {
      results_list_at8[[length(results_list_at8) + 1]] <- result_at8
    } else {
      cat("Insufficient data for structure:", struct, "and diagnosis:", diag, "\n\n")
    }
  }
}

# Combine all AT8 results into a single data frame and write to CSV
all_stats_at8 <- do.call(rbind, results_list_at8)
write.csv(all_stats_at8, "out.dir.kitani/plot/lm_all_stats_GPNMB_at8.csv", row.names = FALSE, quote = F)

