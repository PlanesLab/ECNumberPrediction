# ============================================================
# Reconstruction of results/Case2/Figures/panelA_4plots_final.jpg
# (Coverage/MCC/Precision/Recall across MetaNetX/KEGG/Rhea/ECREACT).
#
# NOTE ON DATA SOURCE: this does NOT read results/Case2/results-DBs/
# panel_data_4plots.xlsx. That file's MetaNetX and KEGG sheets (Panel_A/
# Panel_B) match the original committed image closely, but its Rhea and
# ECREACT sheets (Panel_C/Panel_D) have drifted from it -- e.g. the xlsx
# says SIMMER's ECREACT Precision is 0.92, but the original image clearly
# shows it level with SIMMER's own MCC/Recall (~0.7), no spike. The xlsx
# was evidently edited after panelA_4plots_final.jpg was rendered, and
# Rhea/ECREACT weren't kept in sync.
#
# So every value below was instead extracted directly from the original
# JPG by calibrated pixel measurement: axis gridline pixel-rows detected
# from the y-axis tick labels give the 0.00/1.00 scale, bar-top pixel-rows
# per color per method give the value. Values that came from panels where
# the xlsx was independently confirmed correct (MetaNetX, KEGG) matched
# the xlsx to within ~0.01-0.02 (i.e. pixel-measurement noise); Rhea/
# ECREACT values are pixel-measured only, no xlsx to cross-check against.
# All rounded to 2 decimals, matching the source data's own precision.
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

OUT_ROOT <- "/scratch/jarcagniriv/ECNumberPrediction/results"

method_order <- c("SelenzymeRF", "SIMMER", "Theia", "BEC-Pred")

metric_labels <- c(
  "coverage" = "Coverage",
  "overall_mcc" = "MCC",
  "overall_precision" = "Precision",
  "overall_recall" = "Recall"
)

metric_colors <- c(
  "Coverage"  = "#808080",
  "MCC"       = "#C26683",
  "Precision" = "#B1C266",
  "Recall"    = "#fc8d62"
)

make_panel_a <- function(df, db_name) {
  summary_metrics <- df %>%
    distinct(method, coverage, overall_precision, overall_recall, overall_mcc) %>%
    pivot_longer(
      cols = c(coverage, overall_precision, overall_recall, overall_mcc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      method = factor(method, levels = method_order),
      metric = factor(metric_labels[metric], levels = c("Coverage", "MCC", "Precision", "Recall"))
    )

  ggplot(summary_metrics, aes(x = method, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = metric_colors, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = db_name) +
    theme_minimal(base_size = 28) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 50, r = 20, b = 20, l = 20),
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5)
    )
}

# ============================================================
# Data -- pixel-measured from results/Case2/Figures/panelA_4plots_final.jpg
# (see note above)
# ============================================================

metanetx_df <- data.frame(
  method = method_order,
  coverage = c(0.99, 0.95, 1.00, 1.00),
  overall_mcc = c(0.69, 0.69, 0.83, 0.82),
  overall_precision = c(0.70, 0.68, 0.86, 0.80),
  overall_recall = c(0.70, 0.73, 0.82, 0.82)
)

kegg_df <- data.frame(
  method = method_order,
  coverage = c(1.00, 0.92, 1.00, 1.00),
  overall_mcc = c(0.68, 0.67, 0.79, 0.78),
  overall_precision = c(0.69, 0.68, 0.79, 0.79),
  overall_recall = c(0.68, 0.71, 0.78, 0.78)
)

rhea_df <- data.frame(
  method = method_order,
  coverage = c(1.00, 0.93, 0.99, 1.00),
  overall_mcc = c(0.78, 0.70, 0.95, 0.81),
  overall_precision = c(0.76, 0.74, 0.96, 0.82),
  overall_recall = c(0.79, 0.73, 0.96, 0.82)
)

ecreact_df <- data.frame(
  method = method_order,
  coverage = c(1.00, 0.91, 1.00, 1.00),
  overall_mcc = c(0.74, 0.72, 0.93, 0.94),
  overall_precision = c(0.75, 0.74, 0.94, 0.96),
  overall_recall = c(0.73, 0.72, 0.94, 0.96)
)

panel_A <- make_panel_a(metanetx_df, "MetaNetX")
panel_B <- make_panel_a(kegg_df, "KEGG")
panel_C <- make_panel_a(rhea_df, "Rhea")
panel_D <- make_panel_a(ecreact_df, "ECREACT")

plot_with_legend <- ggplot(
  metanetx_df %>%
    pivot_longer(cols = c(coverage, overall_precision, overall_recall, overall_mcc),
                 names_to = "metric", values_to = "value") %>%
    mutate(
      method = factor(method, levels = method_order),
      metric = factor(metric_labels[metric], levels = c("Coverage", "MCC", "Precision", "Recall"))
    ),
  aes(x = method, y = value, fill = metric)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = metric_colors, name = NULL) +
  theme_minimal(base_size = 28) +
  theme(legend.text = element_text(size = 28), legend.key.size = unit(1.5, "cm"))

legend <- get_legend(plot_with_legend)

combined_plots <- plot_grid(
  panel_A, panel_B,
  panel_C, panel_D,
  labels = c("A", "B", "C", "D"),
  label_size = 40,
  label_fontface = "bold",
  ncol = 2,
  align = "hv"
)

combined_ABCD <- plot_grid(combined_plots, legend, ncol = 2, rel_widths = c(1, 0.15))

ggsave(
  file.path(OUT_ROOT, "Case2/Figures/panelA_4plots_final.jpg"),
  combined_ABCD, width = 22, height = 16, dpi = 300, bg = "white"
)

cat("Saved to Case2/Figures/panelA_4plots_final.jpg\n")
