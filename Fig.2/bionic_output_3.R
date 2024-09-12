rm(list=ls())
options(stringsAsFactors=FALSE)

path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/net_bionic_sig/dim_512_mci.mild"
setwd(path)

library(dplyr)
library(Rtsne)
library(igraph)
library(corrplot)
library(pheatmap)

out.dir <- "/out.dir.kitani"
dir.create(paste0(path, out.dir), recursive = TRUE)

# Load data
bionic <- read.table("dim_512_mci_features.tsv", header = TRUE, row.names = 1, sep = "\t")
bionic_t <- as.data.frame(t(bionic))

# Normalize data
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
data_std <- as.data.frame(lapply(bionic, normalize))

# Distance matrix and adjacency matrix creation
dist_matrix <- proxy::dist(data_std)
set.seed(123)
k <- 25
knn_matrix <- t(apply(dist_matrix, 1, function(x) {
  sort(x, index.return = TRUE)$ix[2:(k + 1)]
}))
adj_matrix <- matrix(0, nrow = nrow(data_std), ncol = nrow(data_std))
for (i in 1:nrow(knn_matrix)) {
  adj_matrix[i, knn_matrix[i, ]] <- 1
}

# Create graph and perform Louvain clustering
graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
communities <- cluster_louvain(graph, resolution = 1)
cluster_result <- membership(communities)

# t-SNE for dimension reduction
set.seed(123)
tsne_result <- Rtsne(data_std, dims = 2)
final_data <- data.frame(tsne_result$Y)
final_data$cluster <- cluster_result
final_data$Gene <- row.names(bionic)
final_data$label <- ifelse(final_data$Gene %in% c("APOE", "MAPT", "APP"), final_data$Gene, NA)

# Save clustering data
write.table(final_data, file = paste0(path, out.dir, "/cliustering_data_k_20.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# Visualization using ggplot2
p <- ggplot(final_data, aes(x = X1, y = X2, color = factor(cluster))) +
  geom_point(alpha = 0.6, size = 3) +
  scale_color_discrete(name = "Cluster") +
  theme_minimal() +
  labs(x = "t-SNE 1", y = "t-SNE 2")

# Save plot as SVG
ggsave(filename = paste0(path, out.dir, "/tsne_plot.svg"), plot = p, device = "svg", width = 12, height = 8)

# Cluster centroids and correlation matrix
genes <- data_std
row.names(genes) <- row.names(bionic)
genes$cluster <- final_data$cluster
centroids <- aggregate(. ~ cluster, data = genes, mean)
cor_matrix <- cor(t(centroids[, -1]))

# Visualize correlation matrix
svg(paste0(path, out.dir, "/corrplot_output.svg"), width = 8, height = 8)
corrplot(cor_matrix, method = "circle")
dev.off()

# Cluster distance matrix and heatmap
dist_matrix_cluster <- dist(centroids)
svg(paste0(path, out.dir, "/pheatmap_output.svg"), width = 8, height = 8)
pheatmap(as.matrix(dist_matrix_cluster))
dev.off()

# Graph plot with edge colors based on distance
adj_matrix <- 1 / as.matrix(dist_matrix_cluster)
graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
E(graph)$weight <- E(graph)$weight / max(E(graph)$weight)
color_func <- colorRampPalette(c("blue", "black", "red"))
E(graph)$color <- color_func(10)[round(E(graph)$weight * 10 + 1)]

# Save graph plot as SVG
svg(paste0(path, out.dir, "/plot_output.svg"), width = 10, height = 10)
plot(graph, 
     edge.width = E(graph)$weight * 10, 
     edge.color = E(graph)$color, 
     vertex.color = "lightblue", 
     vertex.frame.color = NA, 
     vertex.label.family = "Arial")
image.plot(legend.only = TRUE, col = color_func(10), zlim = c(0, 1), 
           legend.args = list(text = "Edge weight", xpd = TRUE))
dev.off()
