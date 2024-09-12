rm(list=ls())  # Clear the workspace
options(stringsAsFactors=FALSE)  # Disable automatic conversion of strings to factors

# Set working directory
path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output"
setwd(path)

# Load necessary libraries
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(openxlsx)
library(gridExtra)

# Output and script directories
out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"

# Load proteomics expression data and convert to dataframe
proteo_exp <- fread("out.dir.kitani/exp_gene_outlierrm_mci.mildAD.txt")
proteo_exp <- as.data.frame(proteo_exp)
row.names(proteo_exp) <- as.character(proteo_exp$V1)
proteo_exp <- proteo_exp[,-1]

# Extract GPNMB expression
gpnmb_expression <- as.numeric(proteo_exp["GPNMB", ])

# Load metadata
metadata <- read.table("out.dir.kitani/meta_outlierrm_with.mildAD.txt")

# Calculate the 35th and 65th percentile thresholds
upper_35_threshold <- quantile(gpnmb_expression, 0.65)
lower_35_threshold <- quantile(gpnmb_expression, 0.35)

gpnmb_expression <- proteo_exp["GPNMB", ]
upper_35_samples <- names(gpnmb_expression[as.numeric(gpnmb_expression) > upper_35_threshold])
lower_35_samples <- names(gpnmb_expression[as.numeric(gpnmb_expression) < lower_35_threshold])

# Extract relevant metadata for upper and lower 35% samples
upper_35_metadata <- metadata[metadata$SampleID %in% upper_35_samples, ]
lower_35_metadata <- metadata[metadata$SampleID %in% lower_35_samples, ]
upper_35_metadata$group <- "Top 35% GPNMB Expression"
lower_35_metadata$group <- "Bottom 35% GPNMB Expression"
combined_data <- rbind(upper_35_metadata, lower_35_metadata)

# Create cross-tables and heatmaps for upper and lower 35% groups
table_data_upper <- table(upper_35_metadata$braaksc, upper_35_metadata$ceradsc)
df_table_upper <- as.data.frame(as.table(table_data_upper))

table_data_lower <- table(lower_35_metadata$braaksc, lower_35_metadata$ceradsc)
df_table_lower <- as.data.frame(as.table(table_data_lower))

# Save heatmaps as SVG
plot_upper <- ggplot(df_table_upper, aes(Var1, Var2)) +
  geom_tile(aes(fill = Freq), color = "white") +
  geom_text(aes(label = sprintf("%d", Freq)), vjust = 1, size = 8) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of braaksc and ceradsc\n(Top 35% GPNMB Expression)",
       x = "Braak SC",
       y = "CERAD SC",
       fill = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

ggsave("out.dir.kitani/plot/GPNMB/heatmap_upper_35.svg", plot = plot_upper, width = 8, height = 6, dpi = 300)

plot_lower <- ggplot(df_table_lower, aes(Var1, Var2)) +
  geom_tile(aes(fill = Freq), color = "white") +
  geom_text(aes(label = sprintf("%d", Freq)), vjust = 1, size = 8) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Heatmap of braaksc and ceradsc\n(Bottom 35% GPNMB Expression)",
       x = "Braak SC",
       y = "CERAD SC",
       fill = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(size = 22),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))

ggsave("out.dir.kitani/plot/GPNMB/heatmap_lower_35.svg", plot = plot_lower, width = 8, height = 6, dpi = 300)

# Test for bias in msex and apoe_genotype
table_msex <- table(combined_data$group, combined_data$msex)
chisq.test(table_msex)
fisher.test(table_msex)

table_apoe <- table(combined_data$group, combined_data$apoe_genotype)
chisq.test(table_apoe)
fisher.test(table_apoe)
print(table_apoe)
prop.table(table_apoe, margin = 1)

# Test for bias in diagnosis
table_diagnosis <- table(combined_data$group, combined_data$diagnosis)
chisq.test(table_diagnosis)
fisher.test(table_diagnosis)

