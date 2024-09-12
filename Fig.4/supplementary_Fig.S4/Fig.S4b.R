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

# Load DEG table
deg.table <- read.csv("out.dir.kitani/community_df.csv", row.names = 1)
names(deg.table) <- c("Symbol","direction")

# Add module scores for each community
unique_directions <- unique(deg.table$direction)
for (i in 1:length(unique_directions)) {
  genes <- deg.table[deg.table$direction == unique_directions[i],]$Symbol
  seurat <- AddModuleScore(seurat, features = list(genes), name = paste0("community_", unique_directions[i]))
}

# Ensure column names are unique and correctly labeled
meta_col_start <- ncol(seurat@meta.data) - length(unique_directions) + 1
meta_col_end <- ncol(seurat@meta.data)
colnames(seurat@meta.data)[meta_col_start:meta_col_end] <- paste0("community_", unique_directions)

# Plot and save DotPlot and FeaturePlot for the communities
svg(file = paste0(out.dir, "/plot_hspa5/seurat_dotplot_output.svg"), width = 8, height = 6)
DotPlot(seurat, features = paste0("community_", unique_directions)) + RotatedAxis() + theme_minimal()
dev.off()

svg(file = paste0(out.dir, "/plot_hspa5/seurat_featureplot_output.svg"), width = 10, height = 5)
FeaturePlot(seurat, features = paste0("community_", unique_directions), ncol = 3) + theme_minimal()
dev.off()

# Save updated Seurat object
saveRDS(seurat, file = "/share1/kitani/data_from_first/GSE174367_scrna.seq/out.dir/rna.seurat.rds")
