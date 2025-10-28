# ================================ #
# Libraries
# ================================ #
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(forcats)
library(stringr)
library(patchwork)
library(cowplot)
library(scales)

# ================================ #
# Data Preparation
# ================================ #
data_eval <- read_csv("/results/Case1/evaluation_summary.csv")

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

# Add class labels
data_eval <- data_eval %>%
  mutate(
    ec_class = as.factor(ec_class),
    ec_class_name = ec_class_names[as.character(ec_class)]
  )

# Keep full dataset for MCC plotting
data_class_mcc <- data_eval

# Method order & metric colors
method_order <- c("E-zyme1", "E-zyme2", "BridgIT", "SelenzymeRF",
                  "SIMMER", "Theia", "BEC-Pred", "CLAIRE")

metric_labels <- c(
  "coverage" = "Coverage",
  "overall_precision" = "Precision",
  "overall_recall" = "Recall",
  "overall_mcc" = "MCC"
)

metric_colors <- c(
  "Coverage" = "#66c2a5",
  "Precision" = "#B1C266",
  "Recall" = "#fc8d62",
  "MCC" = "#C26683"
)

# ================================ #
# Panel A
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
  )

panel_a <- ggplot(summary_metrics, aes(x = method, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = metric_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 28) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.text = element_text(size = 28),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 50, r = 20, b = 20, l = 20)
  )

# ================================ #
# Panel B
# ================================ #
class_mcc_plots <- lapply(method_order, function(m) {
  plot_data <- data_class_mcc %>% filter(method == m)
  ggplot(plot_data, aes(x = ec_class_name, y = class_mcc, fill = ec_class_name)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = class_color_map, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(title = m, x = NULL, y = NULL) +
    theme_minimal(base_size = 20) +
    theme(
      axis.text.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 24, hjust = 0.5, margin = margin(b = 4))
    )
})

legend_plot <- ggplot(data_class_mcc, aes(x = ec_class_name, y = class_mcc, fill = ec_class_name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = class_color_map, name = NULL) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 27),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
legend_b <- get_legend(legend_plot)

panel_b_grid <- wrap_plots(class_mcc_plots, ncol = 4, byrow = TRUE) +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(10, 10, 40, 10))

panel_b_combined <- plot_grid(panel_b_grid, legend_b,
                              ncol = 2, rel_widths = c(5, 1), align = "h")

# ================================ #
# Panel C
# ================================ #
data_topn <- read.csv("/results/Case1/merged_output.csv")

method_order_c <- method_order
method_colors_c <- c(
  "E-zyme1" = "#66C2A5", "E-zyme2" = "#FC8D62", "BridgIT" = "#8DA0CB",
  "SelenzymeRF" = "#E78AC3", "SIMMER" = "#A6D854", "Theia" = "#FFD92F",
  "BEC-Pred" = "#E5C494", "CLAIRE" = "#B3B3B3"
)

extract_subclass <- function(ec_str) {
  if (is.na(ec_str) || !is.character(ec_str)) return(character(0))
  ec_list <- unlist(strsplit(ec_str, "[;|]"))
  subclasses <- str_extract(ec_list, "^\\d+\\.\\d+\\.\\d+")
  unique(na.omit(subclasses))
}

calculate_precision_recall <- function(data, methods, max_top_n = 5) {
  results <- data.frame()
  for (method in methods) {
    for (top_n in 1:max_top_n) {
      precision_list <- c()
      recall_list <- c()
      for (i in seq_len(nrow(data))) {
        true_labels <- extract_subclass(data$EC.Number[i])
        if (length(true_labels) == 0) next
        pred_raw <- data[[method]][i]
        if (is.na(pred_raw) || pred_raw == "") {
          precision_list <- c(precision_list, 0)
          recall_list <- c(recall_list, 0)
          next
        }
        top_preds <- unlist(strsplit(pred_raw, ";"))
        top_preds <- top_preds[1:min(top_n, length(top_preds))]
        all_pred_ecs <- unlist(strsplit(top_preds, "|", fixed = TRUE))
        pred_labels <- extract_subclass(paste(all_pred_ecs, collapse = ";"))
        TP <- sum(pred_labels %in% true_labels)
        FP <- sum(!pred_labels %in% true_labels)
        FN <- sum(!true_labels %in% pred_labels)
        precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
        recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
        precision_list <- c(precision_list, precision)
        recall_list <- c(recall_list, recall)
      }
      results <- rbind(results, data.frame(
        method = method,
        top_n = top_n,
        precision = mean(precision_list, na.rm = TRUE),
        recall = mean(recall_list, na.rm = TRUE)
      ))
    }
  }
  results
}

perf_df <- calculate_precision_recall(data_topn, colnames(data_topn), max_top_n = 5) %>%
  mutate(
    method = recode(method,
                    "E.zyme1" = "E-zyme1",
                    "E.zyme2" = "E-zyme2",
                    "BEC.Pred" = "BEC-Pred")
  ) %>%
  filter(!is.na(method), method %in% method_order_c) %>%
  mutate(method = factor(method, levels = method_order_c))

top_n_labels <- paste0("Top-", 1:5)

panel_c_plot <- function(metric) {
  ggplot(perf_df, aes(x = top_n, y = .data[[metric]], color = method)) +
    geom_line(size = 3) +
    geom_point(size = 5) +
    scale_x_continuous(breaks = 1:5, labels = top_n_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    scale_color_manual(values = method_colors_c, name = NULL) +
    labs(title = stringr::str_to_title(metric), x = NULL, y = NULL) +
    theme_minimal(base_size = 27) +
    theme(
      plot.title = element_text(size = 30, hjust = 0.5, margin = margin(b = 5)),
      legend.position = "none"
    )
}

panel_c <- panel_c_plot("precision")
panel_d <- panel_c_plot("recall")

legend_c <- get_legend(
  ggplot(perf_df, aes(x = top_n, y = precision, color = method)) +
    geom_line(size = 3) +
    geom_point(size = 5) +
    scale_color_manual(values = method_colors_c, name = NULL) +
    theme_minimal(base_size = 27) +
    theme(legend.text = element_text(size = 30))
)

panel_c_combined <- plot_grid(
  plot_grid(panel_c, panel_d, ncol = 2, rel_widths = c(1, 1), align = "h"),
  legend_c, ncol = 2, rel_widths = c(5, 1), align = "h"
)

# ================================ #
# Final combined plot
# ================================ #

# Convert each panel into a patchwork object
panel_a_pw <- wrap_elements(panel_a)
panel_b_pw <- wrap_elements(panel_b_combined)
panel_c_pw <- wrap_elements(panel_c_combined)


# Combine everything first (no labels yet)
final_plot_all <- (panel_a_pw / panel_b_pw / panel_c_pw)

final_plot_labeled <- ggdraw(final_plot_all) +
  draw_plot_label(
    label = c("A", "B", "C"),
    x = c(0, 0, 0),   
    y = c(1, 0.68, 0.35),  
    hjust = -0.5, vjust = 1.2,
    size = 40,         
    fontface = "bold"
  )

ggsave(
  filename = "results/Case1/case1_plot.png",
  plot = final_plot_all,
  width = 20, height = 22, dpi = 300, bg = "white"
)

