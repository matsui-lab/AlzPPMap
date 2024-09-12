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
for (i in 1:length(unique(deg.table$direction))) {
  genes <- deg.table[deg.table$direction == i,]$Symbol
  seurat <- AddModuleScore(seurat, features = list(genes), name = paste0("community_", i))
}

# Rename metadata columns for community scores
names(seurat@meta.data)[19:24] <- paste0("community_", 1:6)

# Plot and save DotPlot and FeaturePlot for the communities
svg("out.dir.kitani/plot_hspa5/seurat_dotplot_output.svg", width = 8, height = 6)
DotPlot(seurat, features = paste0("community_", 1:6)) + RotatedAxis() + theme_minimal()
dev.off()

svg("out.dir.kitani/plot_hspa5/seurat_featureplot_output.svg", width = 10, height = 5)
FeaturePlot(seurat, features = paste0("community_", 1:6), ncol = 3) + theme_minimal()
dev.off()

# Save updated Seurat object
saveRDS(seurat, file = "/share1/kitani/data_from_first/GSE174367_scrna.seq/out.dir/rna.seurat.rds")
