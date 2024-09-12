# Clear environment
rm(list=ls())
options(stringsAsFactors=FALSE)

# Set path and working directory
path <- "/share1/kitani/data_from_first/resilience_ROSMAP_snRNA.seq/microglia"
setwd(path)

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(qs)
library(patchwork)
library(slingshot)
library(BUSpaRse)
library(scales)
library(viridis)
library(tradeSeq)

# Define directories
out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"

# Load Seurat object
seurat <- readRDS("out.dir.kitani/microglia_seurat.rds")

# Run Slingshot for lineage analysis using UMAP embeddings and cluster labels
sds <- slingshot(Embeddings(seurat, "umap"), clusterLabels = seurat$cell_pred)

# Function to generate color palettes for clusters
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# Assign colors to clusters
cell_colors <- cell_pal(seurat$cell_pred, hue_pal())

# Extract dimensional reduction and clustering information
dimred <- seurat@reductions$umap@cell.embeddings
clustering <- seurat$cell_pred

# Run lineage and curve identification using Slingshot
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)

# Define cluster colors
unique_clusters <- unique(clustering)
cluster_colors <- hue_pal()(length(unique_clusters))
names(cluster_colors) <- unique_clusters

# Adjust UMAP plot ranges
x_range <- range(dimred[, 1], na.rm = TRUE)
y_range <- range(dimred[, 2], na.rm = TRUE)
adjusted_x_range <- x_range + 5

# Save UMAP plot as PNG
png("out.dir.kitani/plot/slingshot_plot.png", width = 8, height = 5, units = "in", res = 300)
plot(dimred, col = cell_colors, asp = 1, pch = 20, cex = 0.6, xlab = "UMAP 1", ylab = "UMAP 2", xlim = adjusted_x_range, ylim = y_range)
for (i in seq_along(curve_data)) {
  lines(curve_data[[i]], lwd = ifelse(i == 3, 4, 2), col = "black")
}
legend("topright", legend = names(cluster_colors), col = cluster_colors, pch = 20, title = "Cell Types", cex = 1, pt.cex = 1, bty="n")
dev.off()

# Fetch pseudotime from Slingshot results
pseudotime <- slingPseudotime(sds) %>% as.data.frame()

# Load significant genes
mindy <- read.csv("/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/mindy.app.mapt.csv")
mindy_sig <- mindy[mindy$p.value < 0.05, ]$gene

# Fetch expression data for genes of interest
genes_of_interest <- c(mindy_sig, "GPNMB")
genes_of_interest <- setdiff(genes_of_interest, "GPNMB")
genes_of_interest <- c(genes_of_interest, "GPNMB")
genes_of_interest <- factor(genes_of_interest, levels = genes_of_interest)

gene_expression_data <- FetchData(seurat, vars = genes_of_interest)

# Prepare plot data
df_all <- lapply(genes_of_interest, function(gene) {
  data.frame(
    lineage = rep(c("Lineage 1", "Lineage 2", "Lineage 3"), each = nrow(pseudotime)),
    pseudotime = c(pseudotime$Lineage1, pseudotime$Lineage2, pseudotime$Lineage3),
    gene_expression = gene_expression_data[[gene]],
    gene = gene
  )
}) %>% bind_rows() %>% filter(lineage == "Lineage 3")

# Scale gene expression
df_filtered <- df_all %>%
  group_by(gene) %>%
  mutate(scaled_expression = scale(gene_expression))

# Plot pseudotime vs scaled gene expression
color_palette <- hue_pal()(length(genes_of_interest))
p <- ggplot(df_filtered, aes(x = pseudotime, y = scaled_expression, color = gene)) +
  geom_smooth(method = "gam", se = FALSE, size = 1) +
  scale_color_manual(values = color_palette) +
  facet_wrap(~ gene, scales = "free_y", ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Pseudotime") +
  ylab("Scaled Gene Expression") +
  ggtitle("Pseudotime vs Scaled Gene Expression for Lineage 3")

# Save plot as SVG
svg("./out.dir.kitani/plot/pseudotime_each_gene_expression_plot.svg", width = 4.5, height = 8)
print(p)
dev.off()
