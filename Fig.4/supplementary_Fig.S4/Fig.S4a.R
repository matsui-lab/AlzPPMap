rm(list = ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/net_bionic_sig/dim_512_mci.mild"
setwd(path)

library(dplyr)
library(ggplot2)
library(Seurat)

out.dir <- "/out.dir.kitani"

#############################
#### Single cell Parenchymal
#############################

# Load Seurat object
seurat <- readRDS("/share1/kitani/data_from_first/GSE174367_scrna.seq/out.dir/rna.seurat.rds")
seurat@active.ident <- as.factor(seurat$Cell.Type)

# Create DotPlot and extract data
p <- DotPlot(seurat, features=c("community_1","community_3","community_4","community_5")) + RotatedAxis()
data <- p$data
data$combined_score <- data$avg.exp * data$pct.exp

# Normalize combined scores within each community (feature)
data <- data %>%
  group_by(features.plot) %>%
  mutate(
    combined_score_scaled = (combined_score - min(combined_score)) / (max(combined_score) - min(combined_score))
  ) %>%
  ungroup()

# Plot the normalized combined scores
p <- ggplot(data, aes(x = id, y = combined_score_scaled, fill = combined_score_scaled)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ features.plot) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "ID", y = "Normalized Combined Score") +
  guides(fill = guide_legend(title = NULL)) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))

# Save the plot as SVG
svg("out.dir.kitani/plot_hspa5/module_bar_plot_output.svg", width = 8, height = 4)
print(p)
dev.off()
