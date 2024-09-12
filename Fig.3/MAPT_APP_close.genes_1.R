rm(list=ls())
options(stringsAsFactors=FALSE)

path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/net_bionic_sig/dim_512_mci.mild"
setwd(path)

library(dplyr)
library(ggplot2)
library(ggVennDiagram)

out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"

# Load and normalize data
bionic <- read.table("dim_512_mci_features.tsv", header=TRUE, row.names=1, sep="\t")
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
data_std <- as.data.frame(lapply(bionic, normalize))
row.names(data_std) <- row.names(bionic)

###########################
## Euclidean Distance
###########################

# Calculate Euclidean distance matrix
dist_matrix <- as.matrix(dist(data_std, method = "euclidean"))

# Extract upper triangular part of the matrix
upper_tri <- upper.tri(dist_matrix)
indices <- which(upper_tri, arr.ind = TRUE)

# Create dataframe with distances
result_df <- data.frame(
  Gene1 = rownames(dist_matrix)[indices[, 1]],
  Gene2 = rownames(dist_matrix)[indices[, 2]],
  Distance = dist_matrix[upper_tri]
)
write.table(result_df, paste0(path, out.dir, "/eculidean_all.txt"), quote=FALSE)

# Filter distances for MAPT and APP
mapt_df <- result_df[result_df$Gene1 %in% c("MAPT") | result_df$Gene2 %in% c("MAPT"), ]
app_df <- result_df[result_df$Gene1 %in% c("APP") | result_df$Gene2 %in% c("APP"), ]

# Histograms
hist(mapt_df$Distance, main="Histogram of Gene Distances [MAPT]", 
     xlab="Distance", col="lightblue", breaks=100)
hist(app_df$Distance, main="Histogram of Gene Distances [APP]", 
     xlab="Distance", col="lightblue", breaks=100)

# Filter for 5th percentile
mapt_threshold <- quantile(mapt_df$Distance, 0.05)
filtered_mapt_df <- mapt_df[mapt_df$Distance <= mapt_threshold, ]

app_threshold <- quantile(app_df$Distance, 0.05)
filtered_app_df <- app_df[app_df$Distance <= app_threshold, ]

# Unique genes for MAPT and APP
app_unique <- unique(c(filtered_app_df$Gene1, filtered_app_df$Gene2))
mapt_unique <- unique(c(filtered_mapt_df$Gene1, filtered_mapt_df$Gene2))

# Intersection of MAPT and APP
intersect(app_unique, mapt_unique)

# Venn diagram
lists <- list(APP = app_unique, MAPT = mapt_unique)
ggVennDiagram(lists)

# Remove MAPT and APP from the list and save
filtered_list <- unique(c(app_unique, mapt_unique))
filtered_list <- setdiff(filtered_list, c("MAPT", "APP"))
writeLines(filtered_list, paste0(path, out.dir, "/eculidean_MAPT_APP.txt"))

###########################
## Mediation Analysis
###########################

selected_rows <- data_std[c("MAPT", "APP"), ]
selected_rows <- as.data.frame(t(selected_rows))

# Correlation between MAPT and APP
cor_value <- cor(selected_rows$MAPT, selected_rows$APP)

# Plot correlation
ggplot(selected_rows, aes(x = MAPT, y = APP)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_minimal() + 
  annotate("text", x = Inf, y = Inf, 
           label = sprintf("r = %.2f", cor_value), 
           hjust = 1.1, vjust = 2)
