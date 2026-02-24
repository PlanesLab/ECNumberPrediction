# ============================================================
# Libraries
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(cowplot)
library(scales)

# ============================================================
# Global settings
# ============================================================

method_order <- c(
  "E-zyme1", "E-zyme2", "BridgIT", "SelenzymeRF",
  "SIMMER", "Theia", "BEC-Pred", "CLAIRE"
)

method_colors <- c(
  "E-zyme1"    = "#66C2A5",
  "E-zyme2"    = "#FC8D62",
  "BridgIT"    = "#8DA0CB",
  "SelenzymeRF"= "#E78AC3",
  "SIMMER"     = "#A6D854",
  "Theia"      = "#FFD92F",
  "BEC-Pred"   = "#E5C494",
  "CLAIRE"     = "#B3B3B3"
)

ec_class_names <- c(
  "1" = "Oxidoreductases",
  "2" = "Transferases",
  "3" = "Hydrolases",
  "4" = "Lyases",
  "5" = "Isomerases",
  "6" = "Ligases",
  "7" = "Translocases"
)

class_colors <- c(
  "Oxidoreductases" = "#f3aa7c",
  "Transferases"    = "#a5dbbc",
  "Hydrolases"      = "#8adadf",
  "Lyases"          = "#f3e57c",
  "Isomerases"      = "#E58A98",
  "Ligases"         = "#a5a9db"
)

class_order <- c(
  "Oxidoreductases", "Transferases", "Hydrolases",
  "Lyases", "Isomerases", "Ligases"
)

metric_labels <- c(
  coverage          = "Coverage",
  overall_precision = "Precision",
  overall_recall    = "Recall",
  overall_mcc       = "MCC"
)

metric_colors <- c(
  "Coverage"  = "#66c2a5",
  "Precision" = "#B1C266",
  "Recall"    = "#fc8d62",
  "MCC"       = "#C26683"
)

# ============================================================
# Load data
# ============================================================

data_eval <- read_csv(
  "/Users/josefinaarcagni/Downloads/evaluation_summary_case1.csv",
  show_col_types = FALSE
) %>%
  mutate(
    ec_class_name = ec_class_names[as.character(ec_class)],
    ec_class_name = factor(ec_class_name, levels = class_order),
    method = factor(method, levels = method_order)
  ) %>%
  filter(!is.na(ec_class_name))

data_topn <- read.csv(
  "/Users/josefinaarcagni/Downloads/merged_output_case1.csv",
  stringsAsFactors = FALSE
)

# ============================================================
# Panel A — Overall metrics bar chart
# ============================================================

summary_metrics <- data_eval %>%
  distinct(method, coverage, overall_precision, overall_recall, overall_mcc) %>%
  pivot_longer(
    cols = c(coverage, overall_precision, overall_recall, overall_mcc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = metric_labels[metric])

panel_a <- ggplot(summary_metrics, aes(x = method, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = metric_colors, name = NULL) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 28) +
  theme(
    axis.title   = element_blank(),
    axis.text.x  = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.title = element_blank(),
    legend.text  = element_text(size = 28),
    plot.margin  = margin(t = 50, r = 10, b = 30, l = 10)
  )

# ============================================================
# Panel B — MCC per EC class, one subplot per method
# ============================================================

class_mcc_plots <- lapply(method_order, function(m) {
  ggplot(
    filter(data_eval, method == m),
    aes(x = ec_class_name, y = class_mcc, fill = ec_class_name)
  ) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = class_colors) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = m) +
    theme_minimal(base_size = 20) +
    theme(
      axis.title   = element_blank(),
      axis.text.x  = element_blank(),
      legend.position = "none",
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      plot.margin  = margin(t = 50, r = 10, b = 30, l = 10)
    )
})

# Legend for panel B
legend_b <- get_legend(
  ggplot(data_eval, aes(x = ec_class_name, y = class_mcc, fill = ec_class_name)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = class_colors, name = NULL) +
    theme_minimal(base_size = 24) +
    theme(legend.text = element_text(size = 24))
)

panel_b <- plot_grid(
  plot_grid(plotlist = class_mcc_plots, ncol = 4),
  legend_b,
  ncol = 2,
  rel_widths = c(5, 1)
)

# ============================================================
# Top-N precision / recall calculation
# ============================================================

method_cols <- c(
  "E.zyme1", "E.zyme2", "BridgIT", "SelenzymeRF",
  "SIMMER", "Theia", "BEC.Pred", "CLAIRE"
)

extract_subclass <- function(x) {
  if (is.null(x) || length(x) != 1 || is.na(x) || x == "") return(character(0))
  unique(str_extract(unlist(strsplit(x, "[;|]")), "^\\d+\\.\\d+\\.\\d+"))
}

calculate_precision_recall <- function(data, methods, max_top_n = 5) {
  
  out <- list()
  
  for (m in methods) {
    
    # Find max number of predictions available for this method
    actual_max_n <- max(sapply(data[[m]], function(x) {
      if (is.na(x) || x == "") return(0)
      length(unlist(strsplit(x, ";")))
    }), na.rm = TRUE)
    
    for (n in seq_len(min(actual_max_n, max_top_n))) {
      
      pr <- rr <- numeric()
      
      for (i in seq_len(nrow(data))) {
        
        true <- extract_subclass(data$EC.Number[i])
        if (length(true) == 0) next
        
        all_preds <- unlist(strsplit(data[[m]][i], ";"))
        if (length(all_preds) < n) next
        
        preds <- unique(unlist(lapply(head(all_preds, n), function(pred) {
          unlist(lapply(unlist(strsplit(pred, "\\|")), extract_subclass))
        })))
        preds <- preds[!is.na(preds)]
        
        TP <- sum(preds %in% true)
        FP <- sum(!preds %in% true)
        FN <- sum(!true %in% preds)
        
        pr <- c(pr, ifelse(TP + FP > 0, TP / (TP + FP), 0))
        rr <- c(rr, ifelse(TP + FN > 0, TP / (TP + FN), 0))
      }
      
      out[[length(out) + 1]] <- data.frame(
        method      = m,
        top_n       = n,
        precision   = mean(pr, na.rm = TRUE),
        recall      = mean(rr, na.rm = TRUE)
      )
    }
  }
  
  bind_rows(out) %>%
    mutate(
      method = recode(method,
                      "E.zyme1"  = "E-zyme1",
                      "E.zyme2"  = "E-zyme2",
                      "BEC.Pred" = "BEC-Pred"
      ),
      method = factor(method, levels = method_order)
    )
}

perf_df <- calculate_precision_recall(data_topn, method_cols)

# ============================================================
# Panel C & D — Top-N precision and recall line plots
# ============================================================

plot_topn <- function(metric, title) {
  ggplot(perf_df, aes(x = top_n, y = .data[[metric]], color = method, group = method)) +
    geom_line(size = 3, na.rm = TRUE) +
    geom_point(size = 7, na.rm = TRUE) +
    scale_x_continuous(breaks = 1:5, labels = paste0("Top-", 1:5)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    scale_color_manual(values = method_colors) +
    labs(title = title) +
    theme_minimal(base_size = 26) +
    theme(
      axis.title      = element_blank(),
      legend.position = "none",
      plot.title      = element_text(hjust = 0.5, size = 28, face = "bold"),
      plot.margin     = margin(t = 50, r = 10, b = 30, l = 10)
    )
}

panel_c <- plot_topn("precision", "Precision")
panel_d <- plot_topn("recall", "Recall")

legend_cd <- get_legend(
  ggplot(perf_df, aes(x = top_n, y = precision, color = method)) +
    geom_line(size = 3) +
    scale_color_manual(values = method_colors, name = NULL) +
    theme_minimal(base_size = 26) +
    theme(legend.text = element_text(size = 26))
)

panel_cd <- plot_grid(
  plot_grid(panel_c, panel_d, ncol = 2),
  legend_cd,
  ncol = 2,
  rel_widths = c(5, 1)
)

# ============================================================
# Combine all panels and save
# ============================================================

final_plot <- plot_grid(
  panel_a,
  panel_b,
  panel_cd,
  ncol = 1,
  labels      = c("A", "B", "C"),
  label_size  = 40,
  label_x     = 0.01,
  label_y     = 0.99,
  hjust       = 0,
  vjust       = 1
)

ggsave(
  "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/Case1/fig1_kegg.jpg",
  final_plot,
  width  = 20,
  height = 22,
  dpi    = 300,
  bg     = "white"
)