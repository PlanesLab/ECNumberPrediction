# ============================================================
# Full Case 1 figure (Panels A, B, C) -- kegg500 bootstrap version of
# results/Case1/metrics-case1.R -- same structure, same method order/colors/style,
# fed from the 10-seed kegg500 bootstrap results instead of the single-split
# evaluation, with MajorityVoteCore4 (SelenzymeRF+SIMMER+BEC-Pred+Theia) added.
#
# Panel B uses facet_wrap instead of a manual per-method plot_grid (equivalent
# result); the three panels are combined with cowplot::plot_grid.
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(cowplot)

ROOT <- "/scratch/jarcagniriv/ECNumberPrediction"
OUT_ROOT <- "/scratch/jarcagniriv/ECNumberPrediction/results"

# ============================================================
# Global settings (unchanged from metrics-case1.R, + MajorityVoteCore4)
# ============================================================

method_order <- c(
  "E-zyme1", "E-zyme2", "BridgIT", "SelenzymeRF",
  "SIMMER", "Theia", "BEC-Pred", "MV"
)

method_colors <- c(
  "E-zyme1"    = "#66C2A5",
  "E-zyme2"    = "#FC8D62",
  "BridgIT"    = "#8DA0CB",
  "SelenzymeRF"= "#E78AC3",
  "SIMMER"     = "#A6D854",
  "Theia"      = "#FFD92F",
  "BEC-Pred"   = "#0066FF",
  "MV"         = "#9B30FF"
)

ec_class_names <- c(
  "1" = "Oxidoreductases", "2" = "Transferases", "3" = "Hydrolases",
  "4" = "Lyases", "5" = "Isomerases", "6" = "Ligases", "7" = "Translocases"
)

class_colors <- c(
  "Oxidoreductases" = "#f3aa7c", "Transferases" = "#a5dbbc", "Hydrolases" = "#8adadf",
  "Lyases" = "#f3e57c", "Isomerases" = "#E58A98", "Ligases" = "#a5a9db"
)

class_order <- c("Oxidoreductases", "Transferases", "Hydrolases", "Lyases", "Isomerases", "Ligases")

# ============================================================
# PANEL A -- Overall metrics bar chart, +/- std error bars
# ============================================================

metric_labels <- c(coverage_mean = "Coverage", mcc_mean = "MCC",
                    ppv_mean = "Precision", recall_mean = "Recall")
metric_sd_lookup <- c(Coverage = "coverage_std", MCC = "mcc_std",
                       Precision = "ppv_std", Recall = "recall_std")
metric_colors <- c("Coverage" = "#66c2a5", "MCC" = "#C26683",
                    "Precision" = "#B1C266", "Recall" = "#fc8d62")

data_summary <- read_csv(file.path(OUT_ROOT, "Case1/kegg500_bootstrap_analysis/kegg500_bootstrap_summary.csv"),
                          show_col_types = FALSE) %>%
  mutate(method = recode(method, "MajorityVoteCore4" = "MV")) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  filter(!is.na(method))

summary_means <- data_summary %>%
  select(method, coverage_mean, ppv_mean, recall_mean, mcc_mean) %>%
  pivot_longer(cols = -method, names_to = "metric", values_to = "value") %>%
  mutate(metric = metric_labels[metric])

summary_sds <- data_summary %>%
  select(method, coverage_std, ppv_std, recall_std, mcc_std) %>%
  pivot_longer(cols = -method, names_to = "metric_sd", values_to = "sd") %>%
  mutate(metric = names(metric_sd_lookup)[match(metric_sd, metric_sd_lookup)])

summary_metrics <- summary_means %>%
  left_join(summary_sds %>% select(method, metric, sd), by = c("method", "metric")) %>%
  mutate(metric = factor(metric, levels = c("Coverage", "MCC", "Precision", "Recall")))

panel_a <- ggplot(summary_metrics, aes(x = method, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.85) +
  geom_errorbar(aes(ymin = pmax(0, value - sd), ymax = pmin(1, value + sd)),
                position = position_dodge(width = 0.9), width = 0.25, linewidth = 0.5) +
  scale_fill_manual(values = metric_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  theme_minimal(base_size = 30) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(hjust = 0.5, size = 33),
    axis.text.y = element_text(size = 26),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 38),
    legend.key.size = unit(1.5, "cm"),
    plot.margin = margin(t = 95, r = 10, b = 80, l = 10)
  )

# ============================================================
# PANEL B -- MCC per EC class, faceted per method (facet_wrap instead of
# cowplot::plot_grid over manually-built per-method subplots -- same result),
# +/- std error bars from per-seed per-class scoring (class_metrics_per_seed.py),
# laid out as a 3x3 grid (ncol=3) for the 9 methods.
# ============================================================

data_class <- read_csv(file.path(OUT_ROOT, "Case1/kegg500_bootstrap_analysis/kegg500_bootstrap_class_metrics_summary.csv"),
                        show_col_types = FALSE) %>%
  mutate(
    method = recode(method, "MajorityVoteCore4" = "MV"),
    ec_class_name = factor(ec_class_name, levels = class_order),
    method = factor(method, levels = method_order)
  ) %>%
  filter(!is.na(ec_class_name), !is.na(method))

panel_b <- ggplot(data_class, aes(x = ec_class_name, y = mcc_mean, fill = ec_class_name)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = pmax(0, mcc_mean - mcc_std), ymax = pmin(1, mcc_mean + mcc_std)),
                width = 0.4, linewidth = 0.4) +
  facet_wrap(~ method, ncol = 4) +
  scale_fill_manual(values = class_colors, name = NULL) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(y = "MCC") +
  theme_minimal(base_size = 26) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 34, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 26),
    strip.text = element_text(face = "bold", size = 30),
    legend.position = "right",
    legend.text = element_text(size = 38),
    legend.key.size = unit(1.5, "cm"),
    plot.margin = margin(t = 45, r = 10, b = 80, l = 10)
  )

# ============================================================
# PANEL C -- Top-N precision / recall (pooled across all 10 seeds), faceted
# ============================================================

seed_files <- file.path(ROOT, sprintf("results/Case1/kegg500_bootstrap/seed%d_merged_output.csv", 0:9))
data_topn <- do.call(rbind, lapply(seed_files, read.csv, stringsAsFactors = FALSE))

method_cols <- c("E.zyme1", "E.zyme2", "BridgIT", "SelenzymeRF", "SIMMER", "Theia",
                  "BEC.Pred", "MajorityVoteCore4")

extract_subclass <- function(x) {
  if (is.null(x) || length(x) != 1 || is.na(x) || x == "") return(character(0))
  unique(str_extract(unlist(strsplit(x, "[;|]")), "^\\d+\\.\\d+\\.\\d+"))
}

calculate_precision_recall <- function(data, methods, max_top_n = 5) {
  out <- list()
  for (m in methods) {
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
        TP <- sum(preds %in% true); FP <- sum(!preds %in% true); FN <- sum(!true %in% preds)
        pr <- c(pr, ifelse(TP + FP > 0, TP / (TP + FP), 0))
        rr <- c(rr, ifelse(TP + FN > 0, TP / (TP + FN), 0))
      }
      out[[length(out) + 1]] <- data.frame(method = m, top_n = n, precision = mean(pr, na.rm = TRUE),
                                            recall = mean(rr, na.rm = TRUE), n_support = length(pr))
    }
  }
  bind_rows(out) %>%
    mutate(
      method = recode(method, "E.zyme1" = "E-zyme1", "E.zyme2" = "E-zyme2", "BEC.Pred" = "BEC-Pred",
                       "MajorityVoteCore4" = "MV"),
      method = factor(method, levels = method_order)
    )
}

# A Top-N point is only as trustworthy as how many reactions actually had N ranked candidates
# to average over -- most methods here rarely rank that many (e.g. SIMMER gives >=4 candidates
# for only 6 of ~4900 pooled reactions), so an unfiltered Top-4/Top-5 point can be a wild swing
# (even a "perfect" 1.0) computed from a near-empty bucket rather than a real signal. Drop any
# point with fewer than MIN_SUPPORT reactions behind it instead of silently plotting it.
MIN_SUPPORT <- 30

perf_df_raw <- calculate_precision_recall(data_topn, method_cols)

dropped <- perf_df_raw %>% filter(n_support < MIN_SUPPORT)
if (nrow(dropped) > 0) {
  cat("Dropped low-support Top-N points (n <", MIN_SUPPORT, "):\n")
  print(dropped %>% select(method, top_n, n_support) %>% arrange(method, top_n))
}

perf_df <- perf_df_raw %>%
  filter(n_support >= MIN_SUPPORT) %>%
  pivot_longer(cols = c(precision, recall), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, precision = "Precision", recall = "Recall"))

panel_c <- ggplot(perf_df, aes(x = top_n, y = value, color = method, group = method)) +
  geom_line(linewidth = 2, na.rm = TRUE) +
  geom_point(size = 9, na.rm = TRUE) +
  facet_wrap(~ metric) +
  scale_x_continuous(breaks = 1:5, labels = paste0("Top-", 1:5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_color_manual(values = method_colors, name = NULL) +
  theme_minimal(base_size = 26) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 26),
    strip.text = element_text(face = "bold", size = 36),
    legend.text = element_text(size = 38),
    legend.key.size = unit(1.5, "cm"),
    panel.spacing = unit(3, "cm"),
    plot.margin = margin(t = 35, r = 10, b = 10, l = 10)
  )

# ============================================================
# Combine A + B + C -- each panel's own legend is pulled out and placed in a
# fixed-ratio column (80% plot / 20% legend) so the actual plot areas end up
# exactly the same width across A, B, and C, regardless of how long each
# legend's labels are (e.g. B's "Oxidoreductases" vs A's "MCC").
# ============================================================

legend_a <- get_legend(panel_a)
legend_b <- get_legend(panel_b)
legend_c <- get_legend(panel_c)

row_a <- plot_grid(panel_a + theme(legend.position = "none"), legend_a,
                    ncol = 2, rel_widths = c(0.86, 0.14))
row_b <- plot_grid(panel_b + theme(legend.position = "none"), legend_b,
                    ncol = 2, rel_widths = c(0.8, 0.2))
row_c <- plot_grid(panel_c + theme(legend.position = "none"), legend_c,
                    ncol = 2, rel_widths = c(0.8, 0.2))

final_plot <- plot_grid(
  row_a, row_b, row_c,
  ncol = 1,
  labels = c("A", "B", "C"),
  label_size = 64,
  label_x = 0.01,
  label_y = 1,
  hjust = 0,
  vjust = c(1.05, 1.05, 1.1),
  rel_heights = c(0.36, 0.293, 0.28)
)

ggsave(
  file.path(OUT_ROOT, "Case1/Figures/kegg500_bootstrap_fig1_full.jpg"),
  final_plot, width = 30, height = 32.85, dpi = 300, bg = "white", limitsize = FALSE
)

cat("Saved to Case1/Figures/kegg500_bootstrap_fig1_full.jpg\n")
