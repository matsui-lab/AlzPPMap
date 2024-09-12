rm(list=ls())
options(stringsAsFactors=FALSE)

path <- "/share1/kitani/data_from_first/ROSMAP/rosmap.proteomics/rosmap_output"
setwd(path)

library(dplyr)
library(data.table)
library(reshape2)

out.dir <- "/out.dir.kitani"
script.dir <- "/script.kitani"

# Load proteomics expression data
proteo_exp <- fread("out.dir.kitani/exp_gene_outlierrm_mci.mildAD.txt")
proteo_exp <- as.data.frame(proteo_exp)
row.names(proteo_exp) <- as.character(proteo_exp$V1)
proteo_exp <- proteo_exp[,-1]

#########################
## P-value threshold
#########################

# Compute correlation matrix and p-values
corrMatrix <- rcorr(as.matrix(t(proteo_exp)), type = "pearson")
p_valueMatrix <- corrMatrix$P

# Convert correlation matrix to data frame
df_corr <- melt(corrMatrix$r)
names(df_corr) <- c("from","to","weight")
df_pvalue <- melt(p_valueMatrix)
names(df_pvalue) <- c("from","to","p_value")

# Merge correlation and p-value data frames
df <- merge(df_corr, df_pvalue, by = c("from","to"))
df$p_adjusted <- p.adjust(df$p_value, method = "BH")
df_2 <- df %>% filter(p_adjusted < 0.05, weight > 0) %>%
  mutate(from = ifelse(from < to, from, to), to = ifelse(from < to, to, from)) %>%
  distinct() %>% filter(from != to)

# Min-max normalization
x <- df_2$weight
df_2$weight_01 <- (x - min(x)) / (max(x) - min(x))

#########################
## Set network for BIONIC
#########################

ppi.phy <- fread("~/data.source/string/out.dir.kitani/ppi.phy.exp_gene_750_0_1_norm.txt")

# Prepare for BIONIC
weight_2 <- df_2[,c("from","to","weight","weight_01")]
ppi.phy$type <- "PPI"
weight_2$type <- "Exp"
names(ppi.phy) <- c("name1","name2","weight","norm_weight","type")
names(weight_2) <- c("name1","name2","weight","norm_weight","type")

# Extract common genes
genes_inexpression <- unique(c(weight_2$name1, weight_2$name2))
genes_inppi <- unique(c(ppi.phy$name1, ppi.phy$name2))
genes_forbionic <- unique(intersect(genes_inexpression, genes_inppi))

co_expression_bionic <- weight_2[weight_2$name1 %in% genes_forbionic & weight_2$name2 %in% genes_forbionic,]
ppi_bionic <- ppi.phy[ppi.phy$name1 %in% genes_forbionic & ppi.phy$name2 %in% genes_forbionic,]

# Further filter common genes
genes_forbionic_2 <- unique(intersect(unique(c(co_expression_bionic$name1, co_expression_bionic$name2)),
                                      unique(c(ppi_bionic$name1, ppi_bionic$name2))))

co_expression_bionic_2 <- co_expression_bionic[co_expression_bionic$name1 %in% genes_forbionic_2 &
                                                 co_expression_bionic$name2 %in% genes_forbionic_2,]
ppi_bionic_2 <- ppi_bionic[ppi_bionic$name1 %in% genes_forbionic_2 &
                             ppi_bionic$name2 %in% genes_forbionic_2,]

# Save networks
bionic <- "/net_bionic_sig"
dir.create(paste0(path,out.dir,bionic), recursive = TRUE)
write.table(ppi_bionic_2, file=paste0(path,out.dir,bionic,"/ppi_bionic_mci.mild.txt"),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(co_expression_bionic_2, file=paste0(path,out.dir,bionic,"/co_expression_bionic_mci.mild.txt"),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

