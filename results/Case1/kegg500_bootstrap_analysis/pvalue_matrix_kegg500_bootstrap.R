# ============================================================
# KEGG-500 bootstrap pairwise significance figure -- reconstruction of
# results/Case1/Figures/kegg500_bootstrap_pvalue_matrix.jpg (an upper
# triangular heatmap despite the "matrix" name; its generating script was
# deleted, but the precomputed pairwise FDR table it read from,
# case1_mcc_method_x_seed_withPVALS.csv, survives). Colors sampled directly from the original image's pixels
# (peach ~ FDR 0, blue ~ FDR 1).
# ============================================================

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

ROOT <- "/scratch/jarcagniriv/ECNumberPrediction"
OUT_ROOT <- "/scratch/jarcagniriv/ECNumberPrediction/results"

method_order <- c("E-zyme1", "E-zyme2", "BridgIT", "SelenzymeRF",
                   "SIMMER", "Theia", "BEC-Pred", "MV")

pvals <- read_csv(
  file.path(ROOT, "results/Case1/kegg500_bootstrap_analysis/case1_mcc_method_x_seed_withPVALS.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Method1 = recode(Method1, "MajorityVoteCore4" = "MV"),
    Method2 = recode(Method2, "MajorityVoteCore4" = "MV")
  ) %>%
  rowwise() %>%
  mutate(
    m1_idx = match(Method1, method_order),
    m2_idx = match(Method2, method_order),
    row_method = if (m1_idx < m2_idx) Method1 else Method2,
    col_method = if (m1_idx < m2_idx) Method2 else Method1
  ) %>%
  ungroup() %>%
  transmute(
    Method1 = factor(row_method, levels = method_order),
    Method2 = factor(col_method, levels = method_order),
    pval, FDR
  )

pvals <- pvals %>%
  mutate(
    label = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)),
    stars = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ "ns"
    ),
    cell_label = paste0(label, "\n", stars)
  )

fig <- ggplot(pvals, aes(x = Method2, y = Method1, fill = FDR)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_label), size = 6, lineheight = 0.9) +
  scale_fill_gradient(low = "#F5E9E2", high = "#8B9AD0", limits = c(0, 1),
                       breaks = seq(0, 1, 0.25), name = NULL) +
  scale_x_discrete(position = "top", limits = method_order) +
  scale_y_discrete(limits = rev(method_order)) +
  coord_fixed() +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 0, size = 20),
    axis.text.y = element_text(size = 20),
    panel.grid = element_blank(),
    legend.key.height = unit(1.4, "cm"),
    plot.margin = margin(t = 80, r = 10, b = 10, l = 10)
  )

ggsave(
  file.path(OUT_ROOT, "Case1/Figures/kegg500_bootstrap_pvalue_matrix.jpg"),
  fig, width = 14, height = 12, dpi = 300, bg = "white"
)

cat("Saved to Case1/Figures/kegg500_bootstrap_pvalue_matrix.jpg\n")
