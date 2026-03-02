######################################################################

## Project: BINF6110 Assignment 2
## Script Purpose: To run enrichment analysis
## Date: 25/02/2026
## Name: Gia Ly

######################################################################

BiocManager::install("clusterProfiler")
BiocManager::install("org.Sc.sgd.db")
library(readr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(org.Sc.sgd.db)
library(ggplot2)

outdir <- file.path("results", "r_outputs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

contrast_label <- "mature_vs_early"  # change if needed
res_file <- file.path(outdir, paste0("DE_", contrast_label, "_shrunken.csv"))

res_df <- read_csv(res_file, show_col_types = FALSE)

stopifnot("ORF" %in% colnames(res_df))

res_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(ORF = as.character(ORF))

sig_orf <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(ORF) %>%
  unique()

bg_orf <- res_df %>%
  pull(ORF) %>%
  unique()

cat("Significant genes:", length(sig_orf), "\n")
cat("Background genes:", length(bg_orf), "\n")

ego_bp <- enrichGO(
  gene          = sig_orf,
  universe      = bg_orf,
  OrgDb         = org.Sc.sgd.db,
  keyType       = "ORF",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Save table (even if empty)
write_csv(as.data.frame(ego_bp),
          file.path(outdir, paste0("ORA_GO_BP_", contrast_label, ".csv")))

# Plot only if there are results
go_df <- as.data.frame(ego_bp)

if (nrow(go_df) > 0) {
  
  p_go <- dotplot(ego_bp, showCategory = 15, color = "p.adjust") +
    scale_fill_viridis(
      option = "plasma",
      direction = -1,
      guide = guide_colorbar(reverse = TRUE)
    )
  
  print(p_go)
  
  ggsave(file.path(outdir, paste0("ORA_GO_BP_", contrast_label, ".png")),
         p_go, width = 8, height = 6, dpi = 300)
  
  message("GO ORA complete. Plot and table saved to: ", outdir)
  
} else {
  message("No significant GO BP enrichment found at current cutoffs.")
}

# slayyyy
