rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/resilience_ROSMAP_snRNA.seq/microglia"
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
seurat <- readRDS("out.dir.kitani/microglia_seurat.rds")

# Generate and save combined plot (DimPlot + DotPlot)
dim_plot <- DimPlot(seurat, label = FALSE)
dot_plot <- DotPlot(seurat, features=c("community_1", "community_3", "community_4", "community_5")) + RotatedAxis()
combined_plot <- dim_plot + dot_plot

# Save the combined plot as SVG
svg(file = paste0(out.dir, "/plot/combined_plot.svg"), width = 1100/72, height = 400/72)
print(combined_plot)
dev.off()

# Save the combined plot as high-resolution PNG
width_in_pixels <- (1100/72) * 300
height_in_pixels <- (400/72) * 300
png(file = paste0(out.dir, "/plot/combined_plot.png"), width = width_in_pixels, height = height_in_pixels, res = 300)
print(combined_plot)
dev.off()

# Save the updated Seurat object
saveRDS(seurat, file = paste0(out.dir, "/microglia_seurat.rds"))

# Generate and save cell type combined plot for community_1 (normalized combined score)
p <- DotPlot(seurat, features=c("community_1", "community_2", "community_3", "community_4", "community_5", "community_6", "community_all")) + RotatedAxis()
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
ggsave(file = paste0(out.dir, "/plot/cell_type_combined_plot.svg"), plot = p, device = "svg", width = 8, height = 6)
