rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/early_AD_scRNA.seq"
setwd(path)

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)

# Define output directories
out.dir <- "/out.dir.kitani"
dir.create(paste0(path, out.dir), recursive = TRUE)

# Load the Seurat object
seurat_obj_sub <- readRDS("./out.dir.kitani/microglia.rds")

# Generate and save combined plot (DimPlot + DotPlot)
dim_plot <- DimPlot(seurat_obj_sub, label = FALSE)
dot_plot <- DotPlot(seurat_obj_sub, features=c("community_1","community_2","community_3","community_4",       
                                               "community_5","community_6")) + RotatedAxis()
combined_plot <- dim_plot + dot_plot

# Save the combined plot as SVG and PNG
svg(file = paste0(out.dir, "/microglia_plot/combined_plot.svg"), width = 1100/72, height = 400/72)
print(combined_plot)
dev.off()

png(file = paste0(out.dir, "/microglia_plot/combined_plot.png"), width = 4400, height = 1600, res = 300)
print(combined_plot)
dev.off()

# Generate and save cell type combined plot (community_1 normalized combined score)
p <- DotPlot(seurat_obj_sub, features=c("community_1", "community_2", "community_3", "community_4",       
                                        "community_5", "community_6", "community_all")) + RotatedAxis()
data <- p$data

# Calculate combined score and normalize it
data <- data %>%
  group_by(features.plot) %>%
  mutate(
    combined_score = avg.exp * pct.exp,
    combined_score_scaled = (combined_score - min(combined_score)) / (max(combined_score) - min(combined_score))
  ) %>%
  ungroup()

# Filter for community_1 data
data_filtered <- data %>%
  filter(features.plot == "community_1")

# Plot normalized combined score for community_1
p <- ggplot(data_filtered, aes(x = id, y = combined_score_scaled, fill = combined_score_scaled)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "ID", y = "Normalized Combined Score") +
  guides(fill = guide_legend(title = NULL)) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))

# Save the cell type combined plot as SVG
ggsave(file = paste0(out.dir, "/microglia_plot/cell_type_combined_plot.svg"), plot = p, device = "svg", width = 8, height = 6)

# Calculate cell counts and frequencies
cell_counts <- seurat_obj_sub@meta.data %>%
  group_by(dataset, status, clusters) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(dataset, status) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

cell_counts$status <- factor(cell_counts$status, levels = c("Ctrl", "Abeta", "AbetaTau"))

# Boxplot for cell composition by status
plot <- ggplot(cell_counts, aes(x = status, y = freq, fill = status)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  labs(x = "Cell Subtype", y = "Composition", fill = "Status") +
  facet_wrap(~clusters, scales = "free_y") +
  theme_minimal()

# Save cell composition plot
ggsave("out.dir.kitani/microglia_plot/cell_composition_by_status.png", plot = plot, width = 10, height = 8, dpi = 300)

# Perform ANOVA for each cluster
anova_results <- list()
for(cluster in unique(cell_counts$clusters)) {
  data_cluster <- filter(cell_counts, clusters == cluster)
  
  # ANOVA
  model <- aov(freq ~ status, data = data_cluster)
  anova_summary <- summary(model)
  anova_results[[cluster]] <- anova_summary
}

# TukeyHSD for specific clusters
tukey_results <- data.frame()
clusters_to_test <- c("Microglia_GPNMB_LPL", "Microglia_GPNMB_PLAT")
for(cluster in clusters_to_test) {
  data_cluster <- filter(cell_counts, clusters == cluster)
  model <- aov(freq ~ status, data = data_cluster)
  tukey_test <- TukeyHSD(model)
  
  # Convert Tukey's test results to a data frame
  tukey_df <- as.data.frame(tukey_test$`status`)
  tukey_df$Cluster <- cluster
  tukey_results <- bind_rows(tukey_results, tukey_df)
}

# Convert ANOVA results to a data frame
anova_df <- data.frame()
for(cluster in unique(cell_counts$clusters)) {
  data_cluster <- filter(cell_counts, clusters == cluster)
  model <- aov(freq ~ status, data = data_cluster)
  anova_summary <- summary(model)
  
  # Convert summary to data frame
  anova_result <- as.data.frame(anova_summary[[1]])
  anova_result$Cluster <- cluster
  anova_df <- bind_rows(anova_df, anova_result)
}

# Save ANOVA and Tukey results as CSV files
write.csv(anova_df, "out.dir.kitani/cell_compopsition_anova.csv", quote = F)
write.csv(tukey_results, "out.dir.kitani/cell_compopsition_tukey.csv", quote = F)

print("ANOVA and Tukey HSD tests completed.")
