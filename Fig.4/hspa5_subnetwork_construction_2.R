rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/net_bionic_sig/dim_512_mci.mild"
setwd(path)

library(dplyr)
library(ggplot2)
library(igraph)

out.dir <- "/out.dir.kitani"

# Load distance data
dist_data <- read.table("out.dir.kitani/eculidean_all.txt")

# Extract genes close to HSPA5, APP, and MAPT
df_list <- lapply(c("MAPT", "APP", "HSPA5"), function(gene) {
  subset_data <- dist_data[dist_data$Gene1 == gene | dist_data$Gene2 == gene, ]
  subset_data$target_gene <- gene
  return(subset_data)
})
df_combined <- do.call(rbind, df_list)

# Plot histograms for gene distances
p <- ggplot(df_combined, aes(x=Distance)) +
  geom_histogram(binwidth=0.3, fill="lightblue", color="black") +
  geom_vline(aes(xintercept=quantile(Distance, 0.05)), color="red", linetype="dashed", linewidth=1) +
  facet_wrap(~ target_gene, ncol=3) +
  labs(title="Gene Distance Histograms", x="Distance", y="Count") +
  theme_minimal(base_size=14)
ggsave(filename = "~/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/plot/HSPA5/gene_distance_histograms.svg",
       plot = p, device = "svg", width = 15, height = 6)

# Function to get top 5% genes by distance
get_top_5_percent_by_gene <- function(dist_data, gene_name) {
  related_to_gene <- dist_data[dist_data$Gene1 == gene_name | dist_data$Gene2 == gene_name, ]
  top_5_percent_data <- related_to_gene[order(related_to_gene$Distance),][1:floor(0.05 * nrow(related_to_gene)),]
  return(top_5_percent_data)
}

# Get top 5% closest genes for HSPA5, APP, and MAPT
top_5_percent_HSPA5 <- get_top_5_percent_by_gene(dist_data, "HSPA5")
top_5_percent_APP <- get_top_5_percent_by_gene(dist_data, "APP")
top_5_percent_MAPT <- get_top_5_percent_by_gene(dist_data, "MAPT")

# Find common genes across HSPA5, APP, and MAPT
common_genes <- Reduce(intersect, list(unique(c(top_5_percent_HSPA5$Gene1, top_5_percent_HSPA5$Gene2)),
                                       unique(c(top_5_percent_APP$Gene1, top_5_percent_APP$Gene2)),
                                       unique(c(top_5_percent_MAPT$Gene1, top_5_percent_MAPT$Gene2))))

# Load bionic data and calculate similarity
bionic <- read.table("dim_512_mci_features.tsv", header=T, row.names=1, sep="\t")
similarity_matrix <- as.matrix(bionic) %*% t(bionic)

# Create graph based on top 95 percentile
threshold <- quantile(similarity_matrix, 0.95)
g <- graph_from_adjacency_matrix(similarity_matrix > threshold, mode="undirected", weighted=TRUE)
g <- delete_edges(g, E(g)[ends(g, E(g))[,1] == ends(g, E(g))[,2]])

# Create subnetwork based on common genes
sub_g <- induced_subgraph(g, V(g)[name %in% common_genes])

# Remove isolated nodes and plot subnetwork
sub_g <- delete_vertices(sub_g, V(sub_g)[degree(sub_g) == 0])

# Define node colors
node_colors <- ifelse(V(sub_g)$name %in% c("MAPT", "APP", "HSPA5"), "red", "gray")

# Plot subnetwork with community detection
community <- cluster_leading_eigen(sub_g)
png(paste0(path, "/", out.dir, "/network_cluster_HSPA5_MAPT_APP.png"), width=1600, height=900, res=300)
plot(community, sub_g, vertex.size=5, vertex.label.cex=0.3, layout=layout_with_fr(sub_g),
     vertex.color=node_colors, edge.color=rgb(0,0,0,0.05))
dev.off()

# Save community data
comm_membership <- membership(community)
community_df <- data.frame(Vertex = V(sub_g)$name, Community = comm_membership)
write.csv(community_df, paste0(path, "/", out.dir, "/community_df.csv"), quote = FALSE)

