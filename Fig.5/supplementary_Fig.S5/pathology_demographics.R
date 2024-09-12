rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/elife_2017"
setwd(path)

library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(openxlsx)
library(gridExtra)

# Define output directories
out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"
dir.create(file.path(path, out.dir), recursive = TRUE)
dir.create(file.path(path, script.dir), recursive = TRUE)

# Load data
data <- read.csv("gene_expression_matrix_2016-03-03/columns-samples.csv")
exp <- read.csv("gene_expression_matrix_2016-03-03/fpkm_table_unnormalized.csv", row.names = 1)
gene <- read.csv("gene_expression_matrix_2016-03-03/rows-genes.csv")
data_hip <- data[data$structure_acronym == "HIP", ]
names(exp) <- gsub("X", "", names(exp))

protein_df <- read.csv("ProteinAndPathologyQuantifications.csv")
meta_df <- read.csv("DonorInformation.csv")

# Meta data overview
meta_df$age <- factor(meta_df$age, levels = c("77", "78", "79", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90-94", "95-99", "100+"))

# Define a consistent theme for all plots
custom_theme <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

# Age distribution by dementia status
p1 <- ggplot(meta_df, aes(x = age)) +
  geom_bar() +
  facet_wrap(~ act_demented) +
  labs(title = "Age Distribution by Dementia Status") +
  custom_theme
ggsave(file.path(out.dir, "plot", "age_distribution_by_dementia_status.svg"), plot = p1, width = 8, height = 4)

# Gender distribution by dementia status
p2 <- ggplot(meta_df, aes(x = sex)) +
  geom_bar() +
  facet_wrap(~ act_demented) +
  labs(title = "Gender Distribution by Dementia Status") +
  custom_theme
ggsave(file.path(out.dir, "plot", "gender_distribution_by_dementia_status.svg"), plot = p2, width = 8, height = 4)

# Braak score distribution by dementia status
p3 <- ggplot(meta_df, aes(x = braak)) +
  geom_bar() +
  facet_wrap(~ act_demented) +
  labs(title = "Braak Score Distribution by Dementia Status") +
  custom_theme
ggsave(file.path(out.dir, "plot", "braak_score_distribution_by_dementia_status.svg"), plot = p3, width = 8, height = 4)

# CERAD score distribution by dementia status
p4 <- ggplot(meta_df, aes(x = cerad)) +
  geom_bar() +
  facet_wrap(~ act_demented) +
  labs(title = "CERAD Score Distribution by Dementia Status") +
  custom_theme
ggsave(file.path(out.dir, "plot", "cerad_score_distribution_by_dementia_status.svg"), plot = p4, width = 8, height = 4)

#################
# Statistical tests
#################

# Age: t-test
meta_df$age_num <- as.numeric(as.character(meta_df$age))  # Convert age to numeric
age_no_dementia <- meta_df$age_num[meta_df$act_demented == "No Dementia"]
age_dementia <- meta_df$age_num[meta_df$act_demented == "Dementia"]

if (length(age_no_dementia) > 1 && length(age_dementia) > 1) {
  t_test_result <- t.test(age_no_dementia, age_dementia)
  print(t_test_result)
} else {
  print("Not enough data for t-test on age.")
}

# Gender: Chi-squared test and Fisher's exact test
sex_table <- table(meta_df$sex, meta_df$act_demented)
if (all(sex_table > 0)) {
  chisq_test_sex <- chisq.test(sex_table)
  fisher_test_sex <- fisher.test(sex_table)
  print(chisq_test_sex)
  print(fisher_test_sex)
} else {
  print("Not enough data for chi-squared or Fisher's exact test on sex.")
}

# Braak score: Fisher's exact test
braak_table <- table(meta_df$braak, meta_df$act_demented)
if (all(braak_table > 0)) {
  fisher_test_braak <- fisher.test(braak_table)
  print(fisher_test_braak)
} else {
  print("Not enough data for Fisher's exact test on Braak score.")
}

# CERAD score: Fisher's exact test
cerad_table <- table(meta_df$cerad, meta_df$act_demented)
if (all(cerad_table > 0)) {
  fisher_test_cerad <- fisher.test(cerad_table)
  print(fisher_test_cerad)
} else {
  print("Not enough data for Fisher's exact test on CERAD score.")
}
