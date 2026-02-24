# ================================ #
# Libraries
# ================================ #
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# ================================ #
# Common settings
# ================================ #

method_order <- c("SelenzymeRF", "SIMMER", "Theia", "BEC-Pred", "CLAIRE")

method_colors_c <- c(
  "SelenzymeRF" = "#E78AC3",
  "SIMMER"      = "#A6D854",
  "Theia"       = "#FFD92F",
  "BEC-Pred"    = "#E5C494",
  "CLAIRE"      = "#B3B3B3"
)

levels_order <- c("Class", "Subclass", "Sub-subclass", "Serial Number")

# ================================ #
# Function to build panel
# ================================ #

make_panel <- function(df, title_text, show_legend = FALSE) {
  
  # Fix possible R column name conversions
  names(df) <- gsub("\\.", " ", names(df))
  
  df_long <- df %>%
    pivot_longer(
      cols = -method,
      names_to = "level",
      values_to = "MCC"
    ) %>%
    mutate(
      level = factor(level, levels = levels_order),
      method = factor(method, levels = method_order),
      MCC = ifelse(level == "Serial Number" & !is.na(MCC) & MCC == 0, NA, MCC)
    )
  
  p <- ggplot(df_long, aes(x = level, y = MCC, color = method, group = method)) +
    geom_line(size = 2.5, na.rm = TRUE) +
    geom_point(size = 4, na.rm = TRUE) +
    scale_y_continuous(
      limits = c(0.3, 1),
      breaks = seq(0.3, 1, 0.1),
      expand = c(0, 0)
    ) +
    scale_color_manual(values = method_colors_c, name = NULL) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 22) +
    theme(
      plot.title = element_text(size = 26, hjust = 0.5),
      legend.position = ifelse(show_legend, "right", "none"),
      legend.text = element_text(size = 22)
    )
  
  return(p)
}

# ============================================================
# ======================  DATA SECTION  ======================
# ============================================================

# A) MetaNetX
metanetx_df <- data.frame(
  method = method_order,
  Class = c(0.8072736, 0.8144871, 0.9306479, 0.9259743, 0.6791303),
  Subclass = c(0.7387079, 0.7349119, 0.8851331, 0.8579916, 0.5412765),
  `Sub-subclass` = c(0.6847791, 0.6851262, 0.8334545, 0.8007309, 0.4654971),
  `Serial Number` = c(0.4467011, 0, 0, 0, 0),
  check.names = FALSE
)

# B) Rhea
rhea_df <- data.frame(
  method = method_order,
  Class = c(0.83, 0.80, 0.98, 0.91, 0.67),
  Subclass = c(0.81, 0.743, 0.96, 0.86, 0.54),
  `Sub-subclass` = c(0.78, 0.70, 0.95, 0.81, 0.47),
  `Serial Number` = c(0.45, NA, NA, NA, NA),
  check.names = FALSE
)

# C) KEGG
kegg_df <- data.frame(
  method = method_order,
  Class = c(0.80, 0.79, 0.88, 0.86, 0.63),
  Subclass = c(0.74, 0.72, 0.84, 0.82, 0.54),
  `Sub-subclass` = c(0.68, 0.67, 0.79, 0.78, 0.49),
  `Serial Number` = c(0.31, 0, 0, 0, 0),
  check.names = FALSE
)

# D) ECReact
ecreact_df <- data.frame(
  method = method_order,
  Class = c(0.87, 0.81, 0.98, 0.985, 0.67),
  Subclass = c(0.84, 0.78, 0.95, 0.96, 0.62),
  `Sub-subclass` = c(0.74, 0.72, 0.93, 0.94, 0.55),
  `Serial Number` = c(0.42, NA, NA, NA, NA),
  check.names = FALSE
)

# ============================================================
# Create panels
# ============================================================

panel_A <- make_panel(metanetx_df, "MetaNetX")
panel_B <- make_panel(rhea_df, "Rhea")
panel_C <- make_panel(kegg_df, "KEGG")
panel_D <- make_panel(ecreact_df, "ECREACT", show_legend = TRUE)

# ============================================================
# Extract single legend
# ============================================================

legend <- get_legend(panel_D)

# Remove legend from panel D copy
panel_D <- panel_D + theme(legend.position = "none")

# ============================================================
# Combine 2x2 grid
# ============================================================

grid_only <- plot_grid(
  panel_A, panel_C,
  panel_B, panel_D,
  ncol = 2
)

final_grid <- plot_grid(
  grid_only,
  legend,
  ncol = 2,
  rel_widths = c(1, 0.25)
)

# Add panel labels
final_labeled <- ggdraw(final_grid) +
  draw_plot_label(
    label = c("A", "B", "C", "D"),
    x = c(0, 0.42, 0, 0.42),
    y = c(1, 1, 0.5, 0.5),
    size = 28,
    fontface = "bold"
  )

# ============================================================
# Save
# ============================================================

ggsave(
  "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/Case2/panelB_4datasets.png",
  final_labeled,
  width = 22,
  height = 14,
  dpi = 300,
  bg = "white"
)