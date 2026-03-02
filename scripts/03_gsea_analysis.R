######################################################################

## Project: BINF6110 Assignment 2
## Script Purpose: GSEA GO Analysis
## Date: 25/02/2026
## Name: Gia Ly

######################################################################

library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(readr)
library(viridis)

contrast_label <- "mature_vs_early"
outdir <- "results/r_outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

res_file <- file.path(outdir, paste0("DE_", contrast_label, "_shrunken.csv"))
res_df <- read_csv(res_file, show_col_types = FALSE)

# Ranked list for GSEA (using log2FoldChange)
res_df <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(ORF = as.character(ORF))

gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$ORF
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Total genes in ranked list:", length(gene_list), "\n")

gsea_go <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

gsea_table <- as.data.frame(gsea_go)

write_csv(gsea_table,
          file.path(outdir, paste0("GSEA_GO_BP_", contrast_label, ".csv")))

if (nrow(gsea_table) > 0) {
  
  p_gsea <- dotplot(gsea_go, showCategory = 15, color = "p.adjust") +
    scale_fill_viridis(option = "plasma", direction = -1)
  
  ggsave(file.path(outdir, paste0("GSEA_GO_BP_", contrast_label, ".png")),
         p_gsea, width = 8, height = 6, dpi = 300)
  
  print(p_gsea)
  
  message("GSEA complete. Saved CSV and PNG to: ", outdir)
  
} else {
  message("No significant GSEA terms found at current cutoff.")
}
