# ================================ #
# Libraries
# ================================ #
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(cowplot)

# ================================ #
# Data Preparation
# ================================ #
# Read evaluation summary
data_eval <- read_csv("/Users/josefinaarcagni/Downloads/evaluation_summary_case2_metanetx.csv")

# EC class names & colors
ec_class_names <- c(
  '1' = 'Oxidoreductases',
  '2' = 'Transferases',
  '3' = 'Hydrolases',
  '4' = 'Lyases',
  '5' = 'Isomerases',
  '6' = 'Ligases',
  '7' = 'Translocases'
)

class_color_map <- c(
  'Oxidoreductases' = '#f3aa7c',
  'Transferases'   = '#a5dbbc',
  'Hydrolases'     = '#8adadf',
  'Lyases'         = '#f3e57c',
  'Isomerases'     = '#E58A98',
  'Ligases'        = '#a5a9db',
  'Translocases'   = '#C2C2C2'
)

data_eval <- data_eval %>%
  mutate(
    ec_class = as.factor(ec_class),
    ec_class_name = ec_class_names[as.character(ec_class)]
  )

# Methods & metrics
method_order <- c("SelenzymeRF", "SIMMER", "Theia", "BEC-Pred", "CLAIRE")

metric_labels <- c(
  "coverage" = "Coverage",
  "overall_precision" = "Precision",
  "overall_recall" = "Recall",
  "overall_mcc" = "MCC"
)

metric_colors <- c(
  "Coverage"  = "#66c2a5",
  "Precision" = "#B1C266",
  "Recall"    = "#fc8d62",
  "MCC"       = "#C26683"
)

method_colors_c <- c(
  "SelenzymeRF" = "#E78AC3",
  "SIMMER"      = "#A6D854",
  "Theia"       = "#FFD92F",
  "BEC-Pred"    = "#E5C494",
  "CLAIRE"      = "#B3B3B3"
)

# ================================ #
# Panel A: Summary Metrics Bar Plot
# ================================ #
summary_metrics <- data_eval %>%
  distinct(method, coverage, overall_precision, overall_recall, overall_mcc) %>%
  pivot_longer(
    cols = c(coverage, overall_precision, overall_recall, overall_mcc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    method = factor(method, levels = method_order),
    metric = metric_labels[metric]
  ) %>%
  filter(method %in% method_order)

panel_a <- ggplot(summary_metrics, aes(x = method, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = metric_colors, name = NULL) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 28) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.text = element_text(size = 28),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 50, r = 20, b = 20, l = 20)
  )

# ================================ #
# Panel B: MCC by EC Level
# ================================ #
mcc_df <- data_eval %>%
  select(method, mcc_class, mcc_subclass, mcc_subsubclass, mcc_full) %>%
  group_by(method) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = c(mcc_class, mcc_subclass, mcc_subsubclass, mcc_full),
    names_to = "level",
    values_to = "MCC"
  ) %>%
  mutate(
    level = factor(
      level,
      levels = c("mcc_class", "mcc_subclass", "mcc_subsubclass", "mcc_full"),
      labels = c("Class", "Subclass", "Sub-subclass", "Serial Number")
    ),
    method = factor(method, levels = method_order)
  ) %>%
  group_by(method) %>%
  mutate(
    # If Serial Number MCC = 0, replace with NA to break the line
    MCC = ifelse(level == "Serial Number" & MCC == 0, NA, MCC)
  ) %>%
  ungroup()

panel_b <- ggplot(mcc_df, aes(x = level, y = MCC, color = method, group = method)) +
  geom_line(size = 3, na.rm = TRUE) +
  geom_point(size = 5, na.rm = TRUE) +
  scale_x_discrete(expand = expansion(add = 0.05)) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = method_colors_c, name = NULL) +
  labs(title = "MCC", x = NULL, y = NULL) +
  theme_minimal(base_size = 27) +
  theme(
    plot.title = element_text(size = 30, hjust = 0.5, margin = margin(b = 5)),
    legend.position = "right",
    legend.text = element_text(size = 27)
  )


# ================================ #
# Final combined plot with labels
# ================================ #
final_plot_ab <- plot_grid(
  panel_a, panel_b,
  ncol = 1,
  rel_heights = c(1, 1),
  align = "v"
)

final_plot_labeled <- ggdraw(final_plot_ab) +
  draw_plot_label(
    label = c("A", "B"),
    x = c(0, 0),
    y = c(1, 0.51),
    hjust = -0.5, vjust = 1.2,
    size = 40,
    fontface = "bold"
  )

# ================================ #
# Save
# ================================ #
ggsave(
  filename = "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/Case2/final_metanetx_plot.png",
  plot = final_plot_labeled,
  width = 18, height = 14, dpi = 300, bg = "white"
)
