rm(list = ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/net_bionic_sig/dim_512_mci.mild/out.dir.kitani"
setwd(path)

library(dplyr)
library(ggplot2)
library(rrvgo)
library(RJSONIO)
library(org.Hs.eg.db)

out.dir <- "/out.dir.kitani"

# Load Enrichment table
enrich <- read.csv(file="enrichment_HSPA5_mediator/EnrichmentTable.csv")

# Load GO BP ID
json_file <- "~/data.source/MSigDB/json/c5.go.bp.v2023.1.Hs.json"
json_data <- fromJSON(json_file)

# Convert JSON data to DataFrame
df <- data.frame(Description = names(json_data), exactSource = sapply(json_data, function(x) x$exactSource))

# Merge enrichment data with GO annotations
enrich_go <- merge(enrich[enrich$geneset == "GObp",], df, by="Description")

# Filter significant enrichment results
enrich_go_sig <- enrich_go[enrich_go$qvalue < 0.05 & enrich_go$Count > 4,]

# RRVGO analysis
ecdir <- "enrichment_HSPA5_mediator/rrvgo"
dir.create(ecdir, recursive = T)

simMatrix <- calculateSimMatrix(enrich_go_sig$exactSource, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
scores <- setNames(-log10(enrich_go_sig$qvalue), enrich_go_sig$exactSource)
reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold=0.9, orgdb="org.Hs.eg.db")

# Create heatmap, scatter plot, and wordcloud
heatmapPlot(simMatrix, reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=10)
scatterPlot(simMatrix, reducedTerms)
treemapPlot(reducedTerms)
wordcloudPlot(reducedTerms, min.freq=1, colors="black")
dev.off()

# Merge reduced terms with enrichment data for plotting
df <- merge(enrich_go_sig, reducedTerms, by.x="exactSource", by.y="go")
df$neg_log_qvalue <- -log10(df$qvalue)

# Plot enrichment results
p <- ggplot(df, aes(x = parentTerm, y = neg_log_qvalue, size = Count, color = parentTerm)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +
  geom_text(data = subset(df, neg_log_qvalue > 8.3), aes(label = ID), vjust = -1, size = 1.7, color = "black") +
  labs(title = "Enrichment Analysis", x = "Parent Term", y = "-log10(qvalue)", size = "Gene Count", color = "Parent Term") +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, size = 9), legend.position = "none")
ggsave(filename = "plot_hspa5/enrichment_analysis_all_plot.svg", plot = p, width = 11, height = 9)
