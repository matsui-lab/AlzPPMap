rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output"
setwd(path)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(infotheo)
library(data.table)

out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"

# Load expression data
proteo_exp <- fread("out.dir.kitani/exp_gene_outlierrm_mci.mildAD.txt")
proteo_exp <- as.data.frame(proteo_exp)
row.names(proteo_exp) <- as.character(proteo_exp$V1)
proteo_exp <- proteo_exp[,-1]

# MINDy function for mutual information calculation
mindy <- function(X,Y,Z,B=1000,method="emp",percent=0.35){
  if(is.vector(Z)) Z <- matrix(Z,ncol=1)
  X <- discretize(X)[,1]
  Y <- discretize(Y)[,1]
  n <- length(X)
  m1 <- round(n*percent,0)
  m2 <- round(n*(1-percent),0)
  id1 <- sort.list(Z)[1:m1]
  id2 <- sort.list(Z)[m2:n]
  score <- abs(mutinformation(X[id2],Y[id2],method) - mutinformation(X[id1],Y[id1],method))
  
  null.score <- replicate(B, {
    id1 <- sample(1:n, m1, replace=FALSE)
    id2 <- sample(1:n, m2, replace=FALSE)
    abs(mutinformation(X[id2], Y[id2], method) - mutinformation(X[id1], Y[id1], method))
  })
  
  p.value <-  (sum(score < null.score) + 1) / (B + 1)
  return(list(score=score, perm.score=null.score, p.value=p.value))
}

# Load gene list
gene_list <- readLines("out.dir.kitani/net_bionic_sig/dim_512_mci.mild/out.dir.kitanieculidean_MAPT_APP.txt")
filtered_list <- setdiff(gene_list, c("MAPT", "APP"))

# Define target genes
X_gene <- "APP"
Y_gene <- "MAPT"

# Initialize results dataframe
results_df <- data.frame()
data_std <- as.data.frame(t(proteo_exp))

# Run MINDy for each gene in the list
for (gene in filtered_list) {
  X <- data_std[, X_gene]
  Y <- data_std[, Y_gene]
  Z <- data_std[, gene]
  
  result <- mindy(X, Y, Z)
  
  temp_df <- data.frame(
    gene = gene,
    score = result$score,
    p.value = result$p.value,
    perm.score.mean = mean(result$perm.score),
    perm.score.sd = sd(result$perm.score)
  )
  
  results_df <- rbind(results_df, temp_df)
}

# Save results to CSV
write.csv(results_df, paste0(path, out.dir, "/mindy.app.mapt.csv"), quote=FALSE)

# Calculate -log10(p) for plotting
results_df$neg_log10_p <- -log10(results_df$p.value)

# Plot results
subset_df <- results_df[results_df$p.value < 0.05, ]

p <- ggplot(results_df, aes(x=score, y=neg_log10_p, color=p.value < 0.05)) +
  geom_point(alpha=0.7) + 
  scale_color_manual(values = c("blue", "red")) + 
  labs(title="Score vs -log10(p)", x="Score", y="-log10(p)") +
  theme_minimal() +
  geom_text_repel(data=subset_df, aes(label=gene), box.padding=0.5, 
                  point.padding=0.5, segment.color="grey50", max.overlaps=20) +
  theme(plot.title = element_text(size = 20),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Save plot as SVG
svg(paste0(path, "/out.dir.kitani/plot/mindy_plot_2.svg"), width = 5, height = 7)
print(p)
dev.off()
