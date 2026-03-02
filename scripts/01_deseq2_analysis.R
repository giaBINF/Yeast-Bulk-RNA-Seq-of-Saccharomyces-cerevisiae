######################################################################

## Project: BINF6110 Assignment 2
## Script Purpose: To run DESeq2 on Salmon output data
## Date: 25/02/2026
## Name: Gia Ly

######################################################################

#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("tximport")
#BiocManager::install("txdbmaker")
#BiocManager::install("airway")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("apeglm")
#BiocManager::install("ggrepel")
#BiocManager::install("pheatmap")

library(tidyverse)
library(dplyr)
library(readr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggrepel)
library(pheatmap)
library(pheatmap)

outdir <- file.path("results", "r_outputs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#===========================================================
# 0.0 Loading Metadata
#===========================================================

meta <- read_csv("data/metadata/metadata.csv")

#===========================================================
# 1.0 tximport-ing
#===========================================================

salmon_root <- "results/salmon"

quant_files <- list.files(
  salmon_root,
  pattern = "quant.sf",
  recursive = TRUE,
  full.names = TRUE
)

# Extract sample names from folder names
sample_names <- basename(dirname(quant_files))
sample_names <- sub("_quant$", "", sample_names)

files <- setNames(quant_files, sample_names)

files

# Sanity check
length(files) # should be 9
all(file.exists(files)) # should be true

# Creating tx2gene
gtf_file <- "data/references/S_cerevisiae_R64_annotation.gtf.gz"

txdb <- makeTxDbFromGFF(gtf_file)

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- AnnotationDbi::select(
  txdb,
  keys = k,
  columns = "GENEID",
  keytype = "TXNAME"
)

tx2gene <- tx2gene[!is.na(tx2gene$GENEID), ]

head(tx2gene)
dim(tx2gene)

# tximport-ing
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene
)

# More sanity checks
names(txi)
dim(txi$counts)
head(colnames(txi$counts))

#===========================================================
# 2.0 DESeq2-ing
#===========================================================
meta <- read_csv("data/metadata/metadata.csv", show_col_types = FALSE) %>%
  mutate(
    stage = factor(tolower(stage), levels = c("early","thin","mature")),
    run = sra_accession
  )

meta <- meta %>% filter(run %in% colnames(txi$counts))
meta <- meta[match(colnames(txi$counts), meta$run), ]
stopifnot(all(meta$run == colnames(txi$counts)))

sample_table <- as.data.frame(meta)
rownames(sample_table) <- sample_table$run

# Create dds and running DESeq2

dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_table,
  design = ~ stage
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

resultsNames(dds)

# For mature vs thin, relevel and rerun
dds_thinref <- dds
dds_thinref$stage <- relevel(dds_thinref$stage, ref = "thin")
dds_thinref <- DESeq(dds_thinref)

#===========================================================
# 3.0 Contrast (alpha = 0.05)
#===========================================================
res_thin_vs_early   <- results(dds, contrast = c("stage", "thin", "early"), alpha = 0.05)
res_mature_vs_early <- results(dds, contrast = c("stage", "mature", "early"), alpha = 0.05)

# Mature vs Thin (relevel)
dds_thinref <- dds
dds_thinref$stage <- relevel(dds_thinref$stage, ref = "thin")
dds_thinref <- DESeq(dds_thinref)
res_mature_vs_thin <- results(dds_thinref, contrast = c("stage", "mature", "thin"), alpha = 0.05)

# LFC shrinkage for nicer plots and ranking
resLFC_thin_vs_early   <- lfcShrink(dds, coef = "stage_thin_vs_early", type = "apeglm")
resLFC_mature_vs_early <- lfcShrink(dds, coef = "stage_mature_vs_early", type = "apeglm")
resLFC_mature_vs_thin  <- lfcShrink(dds_thinref, coef = "stage_mature_vs_thin", type = "apeglm")

# Save results tables
write_csv(as.data.frame(resLFC_thin_vs_early)   %>% tibble::rownames_to_column("ORF"),
          file.path(outdir, "DE_thin_vs_early_shrunken.csv"))
write_csv(as.data.frame(resLFC_mature_vs_early) %>% tibble::rownames_to_column("ORF"),
          file.path(outdir, "DE_mature_vs_early_shrunken.csv"))
write_csv(as.data.frame(resLFC_mature_vs_thin)  %>% tibble::rownames_to_column("ORF"),
          file.path(outdir, "DE_mature_vs_thin_shrunken.csv"))

#===========================================================
# 4.0 PCA Plot
#===========================================================
pca_data <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = stage, shape = stage)) +
  geom_point(size = 4) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  labs(color = "Biofilm Stage",
       shape = "Biofilm Stage")

ggsave(file.path(outdir, "PCA_stage.png"), 
       p_pca, width = 7, height = 5, dpi = 300)

#===========================================================
# 6.0 MA Plot
#===========================================================
save_MA <- function(res, filename) {
  png(file.path(outdir, filename), width = 1200, height = 900, res = 150)
  plotMA(res, ylim = c(-10, 10))
  dev.off()
}

save_MA(resLFC_thin_vs_early,   "MA_thin_vs_early.png")
save_MA(resLFC_mature_vs_early, "MA_mature_vs_early.png")
save_MA(resLFC_mature_vs_thin,  "MA_mature_vs_thin.png")

#===========================================================
# 7.0 Volcano Plot
#===========================================================
volcano_plot <- function(res, filename, title) {
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::mutate(group = dplyr::case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
  
  p <- ggplot(df, aes(log2FoldChange, -log10(padj), color = group)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c(
      "Upregulated" = "#D73027",
      "Downregulated" = "#4575B4",
      "Not Significant" = "grey70"
    )) +
    theme_minimal() +
    labs(
      x = "log2 fold change",
      y = "-log10 adjusted p value",
      color = "Group"
    )
  
  ggsave(file.path(outdir, filename), p, width = 7, height = 5, dpi = 300)
}

volcano_plot(resLFC_thin_vs_early,   "Volcano_thin_vs_early.png")
volcano_plot(resLFC_mature_vs_early, "Volcano_mature_vs_early.png")
volcano_plot(resLFC_mature_vs_thin,  "Volcano_mature_vs_thin.png")

#===========================================================
# 8.0 Heatmap of top 20 mature vs early 
#===========================================================
df_me <- as.data.frame(resLFC_mature_vs_early) %>% na.omit()
top20 <- rownames(df_me[order(df_me$padj), ])[1:20]

mat <- assay(vsd)[top20, , drop = FALSE]
ann <- data.frame(stage = colData(vsd)$stage)
rownames(ann) <- colnames(mat)

png(file.path(outdir, "Heatmap_top20_mature_vs_early.png"), width = 1300, height = 1100, res = 150)
pheatmap(mat, scale = "row", annotation_col = ann, show_colnames = FALSE)
dev.off()

message("Done. Outputs saved to: ", outdir)
