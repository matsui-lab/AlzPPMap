rm(list = ls())
options(stringsAsFactors=FALSE)

path <- "~/ROSMAP/rosmap.proteomics/rosmap_output"
setwd(path)

library(dplyr)
library(data.table)
library(ggplot2)

out.dir <- "/out.dir.kitani"
dir.create(paste0(path, out.dir), recursive=TRUE)

# Load normalized expression data
exp <- read.csv("C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv",
                row.names=1)

# Load metadata
meta.tmt <- read.csv("~/ROSMAP/rosmap.metadata/ROSMAP_assay_proteomics_TMTquantitation_metadata.csv", header = TRUE, row.names = 1)
meta.tmt$SampleID <- gsub("_",".",meta.tmt$SampleID)
meta.tmt <- meta.tmt[meta.tmt$SampleID %in% names(exp), ]
meta.tmt$individualID <- gsub("ROSMAP.DLPFC.", "", row.names(meta.tmt))

meta <- read.csv("~/ROSMAP/rosmap.metadata/ROSMAP_clinical.csv", header = TRUE, row.names = 1)
meta.tmt_2 <- meta.tmt[,c("batch","individualID","SampleID")]
meta <- merge(meta, meta.tmt_2, by="individualID")

# PCA
exp_outna <- na.omit(exp)
pca <- prcomp(t(exp_outna), center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x[,1:2])
scores$SampleID <- row.names(scores)
scores <- merge(scores, meta, by="SampleID")

# Outlier removal
scores_2 <- scores[scores$PC1 < -15,]

# Diagnosis classification
scores_2$diagnosis <- ""
scores_2[scores_2$cogdx==1,]$diagnosis <- "NCI"
scores_2[scores_2$cogdx %in% c(2,3),]$diagnosis <- "MCI"
scores_2[scores_2$cogdx %in% c(4,5),]$diagnosis <- "AD"
scores_2[scores_2$cogdx==6,]$diagnosis <- "others"

# Handling age and diagnosis refinement
scores_2[scores_2$age_first_ad_dx=="90+",]$age_first_ad_dx <- 90
scores_2$age_first_ad_dx <- as.numeric(scores_2$age_first_ad_dx)
scores_3 <- scores_2[scores_2$diagnosis %in% c("NCI","MCI","AD"),]
scores_3$diagnosis_2 <- scores_3$diagnosis
scores_3[scores_3$diagnosis_2=="AD" & scores_3$cts_mmse30_lv  >= 21,]$diagnosis_2 <- "mild_AD"

# Extract samples for MCI and mild AD
mci.mildAD <- scores_3[scores_3$diagnosis_2 %in% c("MCI","mild_AD"),]$SampleID
exp_mci.mildAD <- exp[,mci.mildAD]

# Save the expression and metadata
write.table(exp_mci.mildAD, file=paste0(path, out.dir, "/exp_gene_outlierrm_mci.mildAD.txt"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
write.table(scores_3, file=paste0(path, out.dir, "/meta_outlierrm_with.mildAD.txt"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
