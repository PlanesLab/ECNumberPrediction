# ============================================================
# Rhea Splits Comparison figure -- exact layout of the original
# splitting_comparison_seeds.jpg (no generating script survives in the repo,
# reconstructed by inspecting the original image): 3 stacked panels
# (A = MCC, B = Precision, C = Recall), x = method, fill = split
# (Stratified/Time/Scaffold), grouped bars +/- std error bars (Time has a
# single seed, so no error bar), y in [0,1] with 0.25 gridlines, one shared
# legend. MV = MajorityVoteCore4.
#
# Stratified/Scaffold numbers for SelenzymeRF/SIMMER/Theia/BEC-Pred come
# straight from results/Case2/results-splits/seed_summary.csv, unmodified.
# Time (all methods) and MV (all splits) aren't in that file at all, so
# those are recomputed from the raw seed_runs merged predictions using the
# same overall weighted MCC/precision/recall methodology as
# results/Case1/get_metrics.py -- see splitting_comparison_summary.csv.
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(cowplot)

OUT_ROOT <- "/scratch/jarcagniriv/ECNumberPrediction/results"

method_order <- c("SelenzymeRF", "SIMMER", "Theia", "BEC-Pred", "MV")
split_order  <- c("Stratified", "Time", "Scaffold")

split_colors <- c(
  "Stratified" = "#58C0C0",
  "Time"       = "#E0B858",
  "Scaffold"   = "#B888C8"
)

data_summary <- read_csv(
  file.path(OUT_ROOT, "Case2/results-splits/splitting_comparison_summary.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    method = recode(method, "MajorityVoteCore4" = "MV"),
    method = factor(method, levels = method_order),
    split  = factor(split, levels = split_order)
  )

make_panel <- function(mean_col, sd_col, title) {
  d <- data_summary %>%
    transmute(method, split, value = .data[[mean_col]], sd = .data[[sd_col]])

  ggplot(d, aes(x = method, y = value, fill = split)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
    geom_errorbar(
      aes(ymin = pmax(0, value - ifelse(is.na(sd), 0, sd)),
          ymax = pmin(1, value + ifelse(is.na(sd), 0, sd))),
      position = position_dodge(width = 0.9), width = 0.25, linewidth = 0.5
    ) +
    scale_fill_manual(values = split_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 24) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 20),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 26),
      legend.key.size = unit(1.2, "cm"),
      plot.margin = margin(t = 20, r = 10, b = 20, l = 10)
    )
}

panel_a <- make_panel("mcc_mean", "mcc_std", "MCC")
panel_b <- make_panel("ppv_mean", "ppv_std", "Precision")
panel_c <- make_panel("recall_mean", "recall_std", "Recall")

legend <- get_legend(panel_b)

row_a <- plot_grid(panel_a + theme(legend.position = "none"), NULL, ncol = 2, rel_widths = c(0.86, 0.14))
row_b <- plot_grid(panel_b + theme(legend.position = "none"), legend, ncol = 2, rel_widths = c(0.86, 0.14))
row_c <- plot_grid(panel_c + theme(legend.position = "none"), NULL, ncol = 2, rel_widths = c(0.86, 0.14))

final_plot <- plot_grid(
  row_a, row_b, row_c,
  ncol = 1,
  labels = c("A", "B", "C"),
  label_size = 40,
  label_x = 0.01,
  label_y = 1,
  hjust = 0,
  vjust = 1.05
)

ggsave(
  file.path(OUT_ROOT, "Case2/Figures/splitting_comparison_seeds.jpg"),
  final_plot, width = 18, height = 20, dpi = 300, bg = "white"
)

cat("Saved to Case2/Figures/splitting_comparison_seeds.jpg\n")
