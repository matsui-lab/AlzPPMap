rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/resilience_ROSMAP_snRNA.seq/all"
setwd(path)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(patchwork)

# Define output directories
out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"
dir.create(paste0(path, out.dir), recursive = TRUE)
dir.create(paste0(path, script.dir), recursive = TRUE)

# List of regions and corresponding file names
regions <- c("Angular_gyrus", "Entorhinal_cortex", "Hippocampus", "Midtemporal_cortex", "Prefrontal_cortex", "Thalamus")
file_suffix <- "_bar_plot_output.svg"

# Define a function to process each region
process_region <- function(region) {
  # Load Seurat object for the region
  seurat_file <- paste0("out.dir.kitani/", region, ".rds")
  seurat <- readRDS(seurat_file)
  
  # Generate DotPlot and extract data
  p <- DotPlot(seurat, features = c("community_1", "community_3", "community_4", "community_5")) + RotatedAxis()
  data <- p$data
  
  # Calculate combined score and normalized score
  data <- data %>%
    mutate(combined_score = avg.exp * pct.exp) %>%
    group_by(features.plot) %>%
    mutate(combined_score_scaled = (combined_score - min(combined_score)) / (max(combined_score) - min(combined_score))) %>%
    ungroup()
  
  # Plot the normalized combined score
  plot <- ggplot(data, aes(x = id, y = combined_score_scaled, fill = combined_score_scaled)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ features.plot) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(x = "ID", y = "Normalized Combined Score") +
    guides(fill = guide_legend(title = NULL)) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14))
  
  # Save the plot
  svg_filename <- paste0(out.dir, "/plot/", region, file_suffix)
  svg(svg_filename, width = 8, height = 4)
  print(plot)
  dev.off()
}

# Apply the function to all regions
lapply(regions, process_region)

