library(dplyr)
library(readr)

res_df <- read_csv("results/r_outputs/DE_mature_vs_early_shrunken.csv", show_col_types = FALSE)

cand <- res_df %>%
  filter(!is.na(padj), padj < 0.05, !is.na(log2FoldChange)) %>%
  filter(abs(log2FoldChange) >= 1) %>%
  arrange(padj)

top_up <- cand %>% filter(log2FoldChange > 0) %>% slice_head(n = 5)
top_down <- cand %>% filter(log2FoldChange < 0) %>% slice_head(n = 5)

top_up
top_down

# DEG Table

categorize_deg <- function(res_object) {
  
  df <- as.data.frame(res_object) %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05)
  
  total_deg <- nrow(df)
  
  up4  <- sum(df$log2FoldChange >= 2)
  up2  <- sum(df$log2FoldChange >= 1 & df$log2FoldChange < 2)
  dn2  <- sum(df$log2FoldChange <= -1 & df$log2FoldChange > -2)
  dn4  <- sum(df$log2FoldChange <= -2)
  
  return(c(Up4 = up4,
           Up2 = up2,
           Dn2 = dn2,
           Dn4 = dn4,
           Total = total_deg))
}

thin_early  <- categorize_deg(resLFC_thin_vs_early)
mature_thin <- categorize_deg(resLFC_mature_vs_thin)
mature_early<- categorize_deg(resLFC_mature_vs_early)

deg_summary <- data.frame(
  Category = c("Up4", "Up2", "Dn2", "Dn4", "Total DEG"),
  Thin_vs_Early = thin_early,
  Mature_vs_Thin = mature_thin,
  Mature_vs_Early = mature_early
)

deg_summary
