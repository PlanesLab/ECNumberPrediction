library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(cowplot)

# Define method order and colors
method_order <- c("SelenzymeRF", "SIMMER", "Theia", "BEC-Pred", "CLAIRE")

metric_labels <- c(
  "coverage" = "Coverage",
  "overall_mcc" = "MCC",
  "overall_precision" = "Precision",
  "overall_recall" = "Recall"
)

metric_colors <- c(
  "Coverage"  = "#808080",   # Grey for coverage
  "MCC"       = "#C26683",
  "Precision" = "#B1C266",
  "Recall"    = "#fc8d62"
)

make_panel_a <- function(df, db_name) {
  
  # Rename columns to match expected names
  df <- df %>%
    rename(
      method = Model,
      coverage = Coverage,
      overall_mcc = MCC,
      overall_precision = Precision,
      overall_recall = Recall
    )
  
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
  
  panel_a <- ggplot(summary_metrics, aes(x = method, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = metric_colors, name = NULL) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL, title = db_name) +
    theme_minimal(base_size = 28) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",  # Remove individual legends
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 50, r = 20, b = 20, l = 20),
      plot.title = element_text(size = 32, face = "bold", hjust = 0.5)
    )
  
  return(panel_a)
}

# Read Excel file with all sheets
excel_path <- "/Users/josefinaarcagni/Documents/ECMethods/SupplmenetaryStratified/panel_data_4plots.xlsx"

# Read each sheet
df1 <- read_excel(excel_path, sheet = "MetaNetX")
df2 <- read_excel(excel_path, sheet = "KEGG")
df3 <- read_excel(excel_path, sheet = "Rhea")
df4 <- read_excel(excel_path, sheet = "ECREACT")

# Create panels with database names
panel_A <- make_panel_a(df1, "MetaNetX")
panel_B <- make_panel_a(df2, "KEGG")
panel_C <- make_panel_a(df3, "Rhea")
panel_D <- make_panel_a(df4, "ECREACT")

# Extract legend from one of the panels (before removing it)
# Create a temporary plot with legend to extract
df_temp <- df1 %>%
  rename(
    method = Model,
    coverage = Coverage,
    overall_mcc = MCC,
    overall_precision = Precision,
    overall_recall = Recall
  )

summary_metrics_temp <- df_temp %>%
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

plot_with_legend <- ggplot(summary_metrics_temp, aes(x = method, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = metric_colors, name = NULL) +
  theme_minimal(base_size = 28) +
  theme(
    legend.text = element_text(size = 28),
    legend.key.size = unit(1.5, "cm")
  )

# Extract the legend
legend <- get_legend(plot_with_legend)

# Combine panels with single legend
combined_plots <- plot_grid(
  panel_A, panel_B,
  panel_C, panel_D,
  labels = c("A", "B", "C", "D"),
  label_size = 40,
  label_fontface = "bold",
  ncol = 2,
  align = "hv"
)

# Add legend to the right
combined_ABCD <- plot_grid(
  combined_plots,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.15)
)

# Save plot
ggsave(
  "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/Case2/panelA_4plots_final.jpg",
  combined_ABCD,
  width = 22, height = 16, dpi = 300, bg = "white"
)

print("Plot saved successfully!")