# ============================================================
# Libraries
# ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(cowplot)

# ============================================================
# Load and prepare data
# ============================================================
data <- data.frame(
  MCC = c(0.7668409437149193, 0.9537283746326229, 0.7760909322460818,
          0.6993981348825598, 0.8482560789911313, 0.5927238992043349,
          0.35, 0.47, 0.34,
          0.68, 0.78, 0.50,
          0.65, 0.70, 0.49),
  PPV = c(0.7795079807779408, 0.9574582002626094, 0.7800841124554755,
          0.722050221226957, 0.8437786037144461, 0.5794696091685164,
          0.37, 0.51, 0.37,
          0.63, 0.76, 0.52,
          0.69, 0.74, 0.51),
  Recall = c(0.7787934186471663, 0.9547270306258322, 0.7978947368421052,
             0.7093128390596745, 0.863021420518602, 0.6326402016383113,
             0.356, 0.48, 0.34,
             0.69, 0.79, 0.55,
             0.66, 0.73, 0.5081621),
  Split = rep(c("Time", "Stratified", "Scaffold"), 5),
  Method = rep(c("Theia", "BEC-Pred", "CLAIRE", "SelenzymeRF", "SIMMER"), each = 3)
)

# ============================================================
# Method order
# ============================================================
method_order <- c("SelenzymeRF", "SIMMER", "Theia", "BEC-Pred", "CLAIRE")

# Split strategy colors (brighter versions)
split_colors <- c(
  "Stratified" = "#5FC5C1",  # Bright teal
  "Time" = "#E5B85F",        # Bright golden
  "Scaffold" = "#B88DC9"     # Bright mauve
)

# Prepare data
data <- data %>%
  mutate(
    Method = factor(Method, levels = method_order),
    Split = factor(Split, levels = c("Stratified", "Time", "Scaffold"))
  )

# ============================================================
# Theme settings
# ============================================================
theme_custom <- theme_minimal(base_size = 28) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 28),
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold")
  )

# ============================================================
# Create individual plots
# ============================================================

# MCC plot
plot_mcc <- ggplot(data, aes(x = Method, y = MCC, fill = Split)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = split_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "MCC") +
  theme_custom +
  theme(legend.position = "none")

# PPV plot
plot_ppv <- ggplot(data, aes(x = Method, y = PPV, fill = Split)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = split_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "Precision") +
  theme_custom +
  theme(legend.position = "none")

# Recall plot
plot_recall <- ggplot(data, aes(x = Method, y = Recall, fill = Split)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = split_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "Recall") +
  theme_custom +
  theme(legend.position = "none")

# ============================================================
# Extract legend
# ============================================================
legend_plot <- get_legend(
  ggplot(data, aes(x = Method, y = MCC, fill = Split)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = split_colors) +
    theme_minimal(base_size = 26) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 28))
)

# ============================================================
# Combine plots with labels
# ============================================================
final_plot <- plot_grid(
  plot_grid(plot_mcc, plot_ppv, plot_recall, ncol = 1, 
            labels = c("A", "B", "C"),
            label_size = 40,
            label_x = 0.01,
            label_y = 0.99,
            hjust = 0,
            vjust = 1),
  legend_plot,
  ncol = 2,
  rel_widths = c(5, 1)
)

# ============================================================
# Save plot
# ============================================================
ggsave(
  "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/Case2/splitting_comparison.jpg",
  final_plot,
  width = 16,
  height = 18,
  dpi = 300,
  bg = "white"
)

# Display plot
print(final_plot)