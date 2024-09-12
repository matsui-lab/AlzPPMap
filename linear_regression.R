rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "~/elife_2017"
setwd(path)

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)

# Load datasets
protein_df <- read.csv("ProteinAndPathologyQuantifications.csv")
meta_df <- read.csv("DonorInformation.csv")

# Define function to extract summary statistics for IBA1
extract_summary_stats_iba1 <- function(df, demented_status) {
  results <- lapply(split(df, df$structure_acronym), function(sub_df) {
    model <- lm(ihc_tau2_ffpe ~ ihc_a_beta_ffpe * ihc_iba1_ffpe, data=sub_df)
    summary_result <- summary(model)
    coefs <- summary_result$coefficients
    data.frame(
      structure_acronym = unique(sub_df$structure_acronym),
      act_demented = demented_status,
      intercept_estimate = coefs["(Intercept)", "Estimate"],
      intercept_stderror = coefs["(Intercept)", "Std. Error"],
      a_beta_estimate = coefs["ihc_a_beta_ffpe", "Estimate"],
      a_beta_stderror = coefs["ihc_a_beta_ffpe", "Std. Error"],
      iba1_estimate = coefs["ihc_iba1_ffpe", "Estimate"],
      iba1_stderror = coefs["ihc_iba1_ffpe", "Std. Error"],
      interaction_estimate = ifelse("ihc_a_beta_ffpe:ihc_iba1_ffpe" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:ihc_iba1_ffpe", "Estimate"], NA),
      interaction_stderror = ifelse("ihc_a_beta_ffpe:ihc_iba1_ffpe" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:ihc_iba1_ffpe", "Std. Error"], NA),
      r_squared = summary_result$r.squared,
      adj_r_squared = summary_result$adj.r.squared,
      p_value_intercept = coefs["(Intercept)", "Pr(>|t|)"],
      p_value_a_beta = coefs["ihc_a_beta_ffpe", "Pr(>|t|)"],
      p_value_iba1 = coefs["ihc_iba1_ffpe", "Pr(>|t|)"],
      p_value_interaction = ifelse("ihc_a_beta_ffpe:ihc_iba1_ffpe" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:ihc_iba1_ffpe", "Pr(>|t|)"], NA)
    )
  })
  do.call(rbind, results)
}

# Extract summary stats for IBA1 by dementia status
no_df <- protein_df %>% filter(donor_id %in% meta_df %>% filter(act_demented == "No Dementia") %>% pull(donor_id))
yes_df <- protein_df %>% filter(donor_id %in% meta_df %>% filter(act_demented == "Dementia") %>% pull(donor_id))

no_dementia_stats <- extract_summary_stats_iba1(no_df, "No Dementia")
dementia_stats <- extract_summary_stats_iba1(yes_df, "Dementia")

# Combine and export IBA1 summary stats
all_stats_IBA1 <- rbind(no_dementia_stats, dementia_stats)
all_stats_IBA1$factor <- "IBA1"
write.csv(all_stats_IBA1, paste0(path, "/lm_all_stats_IBA1.csv"), row.names=FALSE, quote=FALSE)

# Define function to extract summary statistics for GFAP
extract_summary_stats_gfap <- function(df, demented_status) {
  results <- lapply(split(df, df$structure_acronym), function(sub_df) {
    model <- lm(ihc_tau2_ffpe ~ ihc_a_beta_ffpe * ihc_gfap_ffpe, data=sub_df)
    summary_result <- summary(model)
    coefs <- summary_result$coefficients
    data.frame(
      structure_acronym = unique(sub_df$structure_acronym),
      act_demented = demented_status,
      intercept_estimate = coefs["(Intercept)", "Estimate"],
      intercept_stderror = coefs["(Intercept)", "Std. Error"],
      a_beta_estimate = coefs["ihc_a_beta_ffpe", "Estimate"],
      a_beta_stderror = coefs["ihc_a_beta_ffpe", "Std. Error"],
      gfap_estimate = coefs["ihc_gfap_ffpe", "Estimate"],
      gfap_stderror = coefs["ihc_gfap_ffpe", "Std. Error"],
      interaction_estimate = ifelse("ihc_a_beta_ffpe:ihc_gfap_ffpe" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:ihc_gfap_ffpe", "Estimate"], NA),
      interaction_stderror = ifelse("ihc_a_beta_ffpe:ihc_gfap_ffpe" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:ihc_gfap_ffpe", "Std. Error"], NA),
      r_squared = summary_result$r.squared,
      adj_r_squared = summary_result$adj.r.squared,
      p_value_intercept = coefs["(Intercept)", "Pr(>|t|)"],
      p_value_a_beta = coefs["ihc_a_beta_ffpe", "Pr(>|t|)"],
      p_value_gfap = coefs["ihc_gfap_ffpe", "Pr(>|t|)"],
      p_value_interaction = ifelse("ihc_a_beta_ffpe:ihc_gfap_ffpe" %in% rownames(coefs), coefs["ihc_a_beta_ffpe:ihc_gfap_ffpe", "Pr(>|t|)"], NA)
    )
  })
  do.call(rbind, results)
}

# Extract summary stats for GFAP by dementia status
no_dementia_stats <- extract_summary_stats_gfap(no_df, "No Dementia")
dementia_stats <- extract_summary_stats_gfap(yes_df, "Dementia")

# Combine and export GFAP summary stats
all_stats_GFAP <- rbind(no_dementia_stats, dementia_stats)
all_stats_GFAP$factor <- "GFAP"
write.csv(all_stats_GFAP, paste0(path, "/lm_all_stats_GFAP.csv"), row.names=FALSE, quote=FALSE)

# Define function to extract summary statistics for AT8 (IBA1 and GFAP)
extract_summary_stats_at8 <- function(df, demented_status, marker) {
  results <- lapply(split(df, df$structure_acronym), function(sub_df) {
    model <- lm(ihc_at8_ffpe ~ ihc_a_beta_ffpe * get(marker), data=sub_df)
    summary_result <- summary(model)
    coefs <- summary_result$coefficients
    data.frame(
      structure_acronym = unique(sub_df$structure_acronym),
      act_demented = demented_status,
      intercept_estimate = coefs["(Intercept)", "Estimate"],
      intercept_stderror = coefs["(Intercept)", "Std. Error"],
      a_beta_estimate = coefs["ihc_a_beta_ffpe", "Estimate"],
      a_beta_stderror = coefs["ihc_a_beta_ffpe", "Std. Error"],
      marker_estimate = coefs[paste0("ihc_", marker, "_ffpe"), "Estimate"],
      marker_stderror = coefs[paste0("ihc_", marker, "_ffpe"), "Std. Error"],
      interaction_estimate = ifelse(paste0("ihc_a_beta_ffpe:ihc_", marker, "_ffpe") %in% rownames(coefs), coefs[paste0("ihc_a_beta_ffpe:ihc_", marker, "_ffpe"), "Estimate"], NA),
      interaction_stderror = ifelse(paste0("ihc_a_beta_ffpe:ihc_", marker, "_ffpe") %in% rownames(coefs), coefs[paste0("ihc_a_beta_ffpe:ihc_", marker, "_ffpe"), "Std. Error"], NA),
      r_squared = summary_result$r.squared,
      adj_r_squared = summary_result$adj.r.squared,
      p_value_intercept = coefs["(Intercept)", "Pr(>|t|)"],
      p_value_a_beta = coefs["ihc_a_beta_ffpe", "Pr(>|t|)"],
      p_value_marker = coefs[paste0("ihc_", marker, "_ffpe"), "Pr(>|t|)"],
      p_value_interaction = ifelse(paste0("ihc_a_beta_ffpe:ihc_", marker, "_ffpe") %in% rownames(coefs), coefs[paste0("ihc_a_beta_ffpe:ihc_", marker, "_ffpe"), "Pr(>|t|)"], NA)
    )
  })
  do.call(rbind, results)
}

# Extract summary stats for AT8 (IBA1)
no_dementia_stats_at8 <- extract_summary_stats_at8(no_df, "No Dementia", "iba1")
dementia_stats_at8 <- extract_summary_stats_at8(yes_df, "Dementia", "iba1")
all_stats_IBA1_at8 <- rbind(no_dementia_stats_at8, dementia_stats_at8)
all_stats_IBA1_at8$factor <- "IBA1"
write.csv(all_stats_IBA1_at8, paste0(path, "/lm_all_stats_IBA1_at8.csv"), row.names=FALSE, quote=FALSE)

# Extract summary stats for AT8 (GFAP)
no_dementia_stats_at8 <- extract_summary_stats_at8(no_df, "No Dementia", "gfap")
dementia_stats_at8 <- extract_summary_stats_at8(yes_df, "Dementia", "gfap")
all_stats_GFAP_at8 <- rbind(no_dementia_stats_at8, dementia_stats_at8)
all_stats_GFAP_at8$factor <- "GFAP"
write.csv(all_stats_GFAP_at8, paste0(path, "/lm_all_stats_GFAP_at8.csv"), row.names=FALSE, quote=FALSE)
