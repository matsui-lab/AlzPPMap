rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "~/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/net_bionic_sig/dim_512_mci.mild"
setwd(path)

library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(openxlsx)
library(ggplot2)
library(gridExtra)

out.dir <- "/out.dir.kitani"

# Load community data
deg.table <- read.csv("out.dir.kitani/community_df.csv", row.names = 1)
names(deg.table) <- c("Symbol", "direction")

#####################
## Gene Annotation ##
#####################

# Load gene annotation
patha = "~/data.source"
gfile = paste0(patha, "/genelist_human_2101gene.txt")
gene.table = read.table(gfile, sep="\t", header=T, fill=T)
colnames(gene.table) = c("Gene.ID", "Symbol", "Gene.Name", "Alias", "Ens.ID")

# Merge with community data
deg.table <- left_join(deg.table, gene.table[c("Symbol", "Gene.ID", "Gene.Name", "Alias", "Ens.ID")], by="Symbol")
deg.table <- na.omit(deg.table)

################
## Enrichment ##
################

ecdir <- paste0(path, out.dir, "/enrichment_community_HSPA5")
dir.create(ecdir, recursive = T)

# Define pathways for enrichment analysis
path.gs <- "~/data.source"
path.msig <- paste0(path.gs, "/MSigDB")
gmts <- c("h.all.v2023.1.Hs.entrez", "c2.cp.kegg.v2023.1.Hs.entrez", "c2.cp.reactome.v2023.1.Hs.entrez", "c5.go.bp.v2023.1.Hs.entrez", "c8.all.v2023.1.Hs.entrez")
gmt.names <- c("Hallmark", "KEGG", "Reactome", "GObp", "Cell_Type")

# Perform enrichment analysis for each community
dr.list <- unique(deg.table$direction)
n <- 0
for (dr in dr.list) {
  dr_indx <- deg.table$direction == dr
  genes <- unique(deg.table$Gene.ID[dr_indx])
  if (length(genes) == 0) next
  
  for (k in 1:length(gmts)) {
    gmtfile <- paste0(path.msig, "/", gmts[k], ".gmt")
    gset <- read.gmt(gmtfile)
    egmt <- enricher(genes, TERM2GENE=gset)
    if (is.null(egmt) || sum(egmt@result$p.adjust < 0.05) < 2) next
    
    edox <- setReadable(egmt, 'org.Hs.eg.db', 'ENTREZID')
    cp <- cnetplot(edox, foldChange=genes, cex_label_category=0.8, cex_label_gene=0.8)
    hp <- heatplot(edox, foldChange=genes) + theme(axis.text=element_text(size=6), axis.title=element_text(size=8))
    
    pdf(paste0(ecdir, "/Enrichment_", dr, "_", gmt.names[k], ".pdf"), width=8, height=20)
    grid.arrange(dotplot(egmt, showCategory=20, font.size=8), hp, cp, ncol=1)
    dev.off()
    
    eres <- edox@result
    eres$geneset <- gmt.names[k]
    eres$direction <- dr
    n <- n + 1
    eres.df <- if (n == 1) eres else rbind(eres.df, eres)
  }
}

# Save enrichment results
write.csv(eres.df, file=paste0(ecdir, "/EnrichmentTable.csv"), row.names = F)

#########################
## Specific Enrichment ##
#########################

# Perform specific enrichment for Community 4
j <- 4
dr <- dr.list[j]
dr_indx <- deg.table$direction == dr
genes <- unique(deg.table$Gene.ID[dr_indx])

# Load specific pathway
gmtfile <- paste0(path.msig, "/c5.go.bp.v2023.1.Hs.entrez.gmt")
gset <- read.gmt(gmtfile)
egmt <- enricher(genes, TERM2GENE=gset)

# Save dotplot for Community 4
svg("out.dir.kitani/plot_hspa5/enrich_dotplot_c4.svg", width=8, height=10)
p <- dotplot(egmt, showCategory=15, font.size=14) + ggtitle("Community 4") +
  theme(plot.title = element_text(size = 24), axis.text.y = element_text(size = 12))
print(p)
dev.off()
