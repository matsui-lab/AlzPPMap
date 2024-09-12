rm(list=ls())
options(stringsAsFactors=FALSE)

path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output/out.dir.kitani/net_bionic_sig/dim_512_mci.mild/out.dir.kitani"
setwd(path)

library(dplyr)
library(ggplot2)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(Seurat)

###############
## Load Data ##
###############

# Load clustering data
deg.table <- read.table("cliustering_data_k_20.txt", header = TRUE)
deg.table$direction <- deg.table$cluster
names(deg.table)[4] <- "Symbol"
deg.table <- deg.table[, c(4, 3, 6)]

# Load gene annotation
patha = "~/data.source"
gene.file = paste0(patha, "/genelist_human_2101gene.txt")
gene.table = read.table(gene.file, sep="\t", quote = "", header = TRUE, fill = TRUE)
colnames(gene.table) = c("Gene.ID", "Symbol", "Gene.Name", "Alias", "Ens.ID")
deg.table = left_join(deg.table, gene.table[c("Symbol", "Gene.ID", "Gene.Name", "Alias", "Ens.ID")], by="Symbol")
deg.table = na.omit(deg.table)

##################
## Enrichment Analysis ##
##################

# Create directory for enrichment results
ecdir = paste0(path, "/enrichment_k_20")
dir.create(ecdir, recursive = TRUE)

# Enrichment analysis using MSigDB gene sets
path.msig = paste0("~/data.source/MSigDB")
gmts = c("h.all.v2023.1.Hs.entrez", "c2.cp.kegg.v2023.1.Hs.entrez", "c2.cp.reactome.v2023.1.Hs.entrez", "c5.go.bp.v2023.1.Hs.entrez")
gmt.names = c("Hallmark", "KEGG", "Reactome", "GObp")

dr.list <- unique(deg.table$direction)
n = 0
for (dr in dr.list) {
  genes <- unique(deg.table$Gene.ID[deg.table$direction == dr])
  if (length(genes) == 0) next
  
  for (k in 1:length(gmts)) {
    gmtfile <- paste0(path.msig, "/", gmts[k], ".gmt")
    gset <- read.gmt(gmtfile)
    egmt <- enricher(genes, TERM2GENE=gset)
    if (is.null(egmt) || sum(egmt@result$p.adjust < 0.05) < 2) next
    
    # Save enrichment plots
    pdf(paste0(ecdir, "/Enrichment_", dr, "_", gmt.names[k], ".pdf"), width=8, height=20)
    grid.arrange(dotplot(egmt, showCategory=20), heatplot(egmt, foldChange=genes), cnetplot(setReadable(egmt, 'org.Hs.eg.db', 'ENTREZID')), ncol=1)
    dev.off()
    
    # Save enrichment results
    eres = setReadable(egmt, 'org.Hs.eg.db', 'ENTREZID')@result
    eres$geneset <- gmt.names[k]
    eres$direction <- dr
    n = n + 1
    if (n == 1) {
      eres.df = eres
    } else {
      eres.df = rbind(eres.df, eres)
    }
  }
}
write.csv(eres.df, file=paste0(ecdir, "/EnrichmentTable.csv"), row.names=FALSE)

#####################
## Seurat Analysis ##
#####################

# Load Seurat object (Parenchymal)
seurat <- readRDS("~/GSE174367_scrna.seq/out.dir/rna.seurat.rds")
seurat@active.ident <- as.factor(seurat$Cell.Type)

# Add module scores
module <- read.table("cliustering_data_k_20.txt", header = TRUE)
for (i in 1:length(unique(module$cluster))) {
  genes <- module[module$cluster == i, ]$Gene
  seurat <- AddModuleScore(seurat, features = list(genes), name = paste0("module_", i))
}
names(seurat@meta.data)[19:46] <- paste0("module_", 1:28)

# Save Seurat plots
svg("plot_seurat/dimplot_output.svg", width = 10, height = 5)
DimPlot(seurat, label = TRUE) + DimPlot(seurat, group.by = "Diagnosis") %>% print
dev.off()

svg("plot_seurat/dotplot_output.svg", width = 10, height = 5)
DotPlot(seurat, features = paste0("module_", 1:28)) + RotatedAxis() %>% print
dev.off()

