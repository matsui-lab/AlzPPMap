rm(list=ls())
options(stringsAsFactors=FALSE)

# Set working directory
path <- "~/ROSMAP/rosmap.proteomics/rosmap_output"
setwd(path)

library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(data.table)
library(openxlsx)

out.dir <- "/out.dir.kitani"

# Load MINDy results and filter significant genes
mindy <- read.csv("out.dir.kitani/mindy.app.mapt.csv", row.names = 1)
mindy.sig <- mindy[mindy$p.value <= 0.05,]$gene

# Load tauome data and extract gene list
tauome <- read.xlsx("~/ptau_proteome_PMID.32812023/sheet4_intersect_ptau_ptorein.xlsx", header = FALSE)
names(tauome) <- tauome[3,]
tauome <- as.data.frame(tauome[-c(1,2,3),])
tauome.gene <- unique(tauome[,4])

# Load amyloid data and extract gene list
amyloid <- read.xlsx("~/amyloidome_PMID.30664241/app_interactome.xlsx")
amyloid.gene <- unique(amyloid$AARS2)

# Find intersection between MINDy significant genes and tauome/amyloid gene lists
intersect_tau <- intersect(mindy.sig, tauome.gene)
intersect_amyloid <- intersect(mindy.sig, amyloid.gene)

# Combine lists for Venn diagram
lists <- list(Mindy = mindy.sig, Amyloid = amyloid.gene, Tau = tauome.gene)

# Plot and save Venn diagram
venn <- ggVennDiagram(lists)
ggsave(filename = "out.dir.kitani/plot/interactome_venn_diagram.svg", plot = venn, device = "svg", width = 8, height = 6)
