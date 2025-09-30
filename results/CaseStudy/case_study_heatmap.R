# === Parameters ===
include_majority_vote <- TRUE   # set to FALSE to exclude TOP1/TOP5 majority vote from plots

# === Load libraries ===
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(data.tree)
library(patchwork)
library(RColorBrewer)
library(tibble)
library(purrr)
library(ggpubr)
library(readxl)
library(readr)

# READ DATAFRAME HERE
merged_output <- read_csv("/Users/josefinaarcagni/Downloads/merged_output_drugs2.csv")

# === Functions ===
collapse_to_third_level <- function(predictions) {
  if (is.na(predictions) || predictions == "" || tolower(predictions) == "nan")
    return(NA)
  preds <- str_split(predictions, ";")[[1]]
  preds <- lapply(preds, function(p) {
    sub_preds <- str_split(p, "\\|")[[1]]
    sub_preds <- str_trim(sub_preds)
    sub_preds <- unique(sub_preds[sub_preds != ""])
    sub_preds <- unique(sapply(sub_preds, function(sp) {
      paste(str_split(sp, "\\.")[[1]][1:3], collapse=".")
    }))
    paste(sub_preds, collapse="|")
  })
  paste(preds, collapse=";")
}

extract_top_ec_levels <- function(predictions, top_n = 5) {
  if (is.na(predictions) || predictions == "" || tolower(predictions) == "nan")
    return(rep(NA, top_n))
  preds <- str_split(predictions, ";")[[1]]
  preds <- unlist(lapply(preds, function(p) str_split(p, "\\|")[[1]]))
  preds <- unique(str_trim(preds))
  preds[1:min(length(preds), top_n)]
}

classify_hit <- function(preds, true_ec) {
  if (all(is.na(preds))) return("No prediction")
  true_prefix <- paste(str_split(true_ec, "\\.")[[1]][1:3], collapse=".")
  for (p in preds) {
    parts <- str_split(p, "\\.")[[1]]
    if (length(parts) >= 3 && paste(parts[1:3], collapse=".") == true_prefix) {
      return(ifelse(p == preds[1], "Top 1", "Top 5"))
    }
  }
  "No hit"
}

score_prediction <- function(preds, true_ec) {
  if (all(is.na(preds))) return(0)
  true_prefix <- paste(str_split(true_ec, "\\.")[[1]][1:3], collapse=".")
  for (i in seq_along(preds)) {
    parts <- str_split(preds[[i]], "\\.")[[1]]
    if (length(parts) >= 3 && paste(parts[1:3], collapse=".") == true_prefix)
      return(6 - i)
  }
  0
}

# === Prepare data ===
ec_methods <- c("E-zyme1","E-zyme2", "BridgIT", "SelenzymeRF","SIMMER", "Theia" ,"BEC-Pred","CLAIRE")
df <- merged_output %>% rename(Theia = Theia_ECREACT)
df <- df %>% rowwise() %>% mutate(across(all_of(ec_methods), ~ collapse_to_third_level(.x))) %>% ungroup()

# === Heatmap data ===
heatmap_data <- df %>%
  mutate(drug_EC = str_trim(drug_EC)) %>%
  select(drug, drug_EC, all_of(ec_methods)) %>%
  pivot_longer(all_of(ec_methods), names_to="method", values_to="predictions") %>%
  rowwise() %>%
  mutate(
    top_preds = list(extract_top_ec_levels(predictions, 5)),
    hit_type = classify_hit(top_preds, drug_EC),
    top1_score = as.integer(hit_type == "Top 1"),
    weighted_score = score_prediction(top_preds, drug_EC)
  ) %>%
  ungroup() %>%
  mutate(
    method = factor(method, levels = ec_methods),
    hit_type = factor(hit_type, levels = c("Top 1", "Top 5", "No hit", "No prediction"))
  )

# === Order drugs ===
tree_data_ordered <- df %>%
  distinct(drug, drug_EC) %>%
  separate(drug_EC, into=c("class","subclass","subsubclass","rest"), sep="\\.", fill="right") %>%
  mutate(across(c(class, subclass, subsubclass), as.character)) %>%
  mutate(
    class_num = as.integer(class),
    subclass_num = as.numeric(subclass),
    subsubclass_num = as.numeric(subsubclass),
    class = paste0( class),
    subclass = paste0(class, ".", subclass),
    subsubclass = paste0(subclass, ".", subsubclass)
  ) %>% arrange(class_num, subclass_num, subsubclass_num, drug)

ordered_drugs <- tree_data_ordered$drug
heatmap_data <- heatmap_data %>% mutate(drug=factor(drug, levels=rev(ordered_drugs)))

# === Annotation bar ===
custom_class_colors <- c("1"="#ef8d50","2"="#47b275","3"="#4ac6cd","4"="#ecd632")
get_shade <- function(col, n) tail(colorRampPalette(c("white", col))(n+1), -1)

subsubclass_colors <- tree_data_ordered %>%
  distinct(class, subsubclass) %>%
  group_by(class) %>%
  mutate(shade = get_shade(custom_class_colors[class[1]], n())) %>%
  ungroup() %>%
  select(subsubclass, shade) %>%
  deframe()

annotation_bar <- tree_data_ordered %>%
  select(drug, subsubclass) %>%
  mutate(drug=factor(drug, levels=rev(ordered_drugs)), color=subsubclass_colors[subsubclass])

label_positions <- annotation_bar %>%
  group_by(subsubclass) %>%
  summarise(y=mean(as.numeric(drug)), .groups="drop") %>%
  left_join(annotation_bar %>% distinct(subsubclass, color), by="subsubclass")

annotation_plot <- ggplot(annotation_bar, aes(y=drug, x=0.95)) +
  geom_tile(aes(fill=color), width=0.003) +
  geom_text(data=label_positions, aes(y=y, x=0.944, label=subsubclass), color="black", hjust=0.1, size=3.3) +
  scale_fill_identity() +
  theme_void() +
  theme(plot.margin=margin(5,5,5,10))

# === Heatmap ===
heatmap <- ggplot(heatmap_data, aes(x=method, y=drug, fill=hit_type)) +
  geom_tile(color="white", size=0.7) +
  scale_fill_manual(values=c("Top 1"="#2ecc71","Top 5"="#3498db","No hit"="#e74c3c","No prediction"="#cccccc")) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position="right",
        legend.text = element_text(size = 15),
        plot.margin=margin(10,10,10,0)) +
  labs(fill=NULL)

# === Grouped Barplots (per method) ===
method_hits5 <- heatmap_data %>% filter(hit_type %in% c("Top 1", "Top 5")) %>%
  left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
  count(method, class) %>% mutate(HitType="Top 5")

method_hits1 <- heatmap_data %>% filter(hit_type=="Top 1") %>%
  left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
  count(method, class) %>% mutate(HitType="Top 1")

top_hits <- bind_rows(method_hits1, method_hits5)

# === Majority vote (optional) ===
if (include_majority_vote) {
  majority_methods <- c("Theia", "BEC-Pred", "SIMMER", "SelenzymeRF")
  
  majority_preds <- df %>%
    select(drug, drug_EC, all_of(majority_methods)) %>%
    rowwise() %>%
    mutate(
      true_prefix = paste(str_split(drug_EC, "\\.")[[1]][1:3], collapse="."),
      
      pred_top1 = {
        weighted_votes <- list()
        for (m in majority_methods) {
          pred_str <- get(m)
          if (!is.na(pred_str)) {
            preds <- str_split(pred_str, ";")[[1]][1]
            if (length(preds) > 0) {
              pipes <- str_split(preds, "\\|")[[1]]
              pipes <- str_trim(pipes)
              pipes <- pipes[pipes != ""]
              for (p in pipes) {
                subsub <- paste(str_split(p, "\\.")[[1]][1:3], collapse=".")
                weighted_votes[[subsub]] <- (weighted_votes[[subsub]] %||% 0) + 5
              }
            }
          }
        }
        if (length(weighted_votes) == 0) NA else names(sort(unlist(weighted_votes), decreasing=TRUE))[1]
      },
      
      pred_top5 = {
        weighted_votes <- list()
        for (m in majority_methods) {
          pred_str <- get(m)
          if (!is.na(pred_str)) {
            preds <- str_split(pred_str, ";")[[1]]
            preds <- preds[1:min(5, length(preds))]
            for (i in seq_along(preds)) {
              pipes <- str_split(preds[[i]], "\\|")[[1]]
              pipes <- str_trim(pipes)
              pipes <- pipes[pipes != ""]
              weight <- max(6 - i, 1)
              for (p in pipes) {
                subsub <- paste(str_split(p, "\\.")[[1]][1:3], collapse=".")
                weighted_votes[[subsub]] <- (weighted_votes[[subsub]] %||% 0) + weight
              }
            }
          }
        }
        if (length(weighted_votes) == 0) NA else names(sort(unlist(weighted_votes), decreasing=TRUE))[1]
      },
      
      top1_hit_type = case_when(
        is.na(pred_top1) ~ "No prediction",
        pred_top1 == true_prefix ~ "Top 1",
        TRUE ~ "No hit"
      ),
      
      top5_hit_type = case_when(
        is.na(pred_top5) ~ "No prediction",
        pred_top5 == true_prefix ~ "Top 5",
        TRUE ~ "No hit"
      )
    ) %>% ungroup()
  
  # summaries for plotting
  top1_hits_only <- heatmap_data %>%
    filter(hit_type == "Top 1") %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    count(method, class) %>%
    mutate(HitType = "Top 1")
  
  majority_only_top1 <- majority_preds %>%
    filter(top1_hit_type == "Top 1") %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    count(method = "TOP1", class) %>%
    mutate(HitType = "TOP1")
  
  majority_only_top5 <- majority_preds %>%
    filter(top5_hit_type %in% c("Top 1", "Top 5")) %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    count(method = "TOP5", class) %>%
    mutate(HitType = "TOP5")
  
  top1_plus_majority <- bind_rows(top1_hits_only, majority_only_top1, majority_only_top5)
  
} else {
  # if not including majority vote
  top1_plus_majority <- heatmap_data %>%
    filter(hit_type == "Top 1") %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    count(method, class) %>%
    mutate(HitType = "Top 1")
}

# set factor levels depending on inclusion
if (include_majority_vote) {
  top1_plus_majority <- top1_plus_majority %>%
    mutate(method = factor(method, levels = c(ec_methods, "TOP1", "TOP5")))
} else {
  top1_plus_majority <- top1_plus_majority %>%
    mutate(method = factor(method, levels = ec_methods))
}

# === Barplot for part B ===
ec_labels <- c(
  "1" = "Oxidoreductases",
  "2" = "Transferases",
  "3" = "Hydrolases",
  "4" = "Lyases"
)

top1_barplot_with_majority <- ggplot(top1_plus_majority, aes(y = n, x = method, fill = class)) +
  geom_hline(yintercept = seq(0, 25, by = 5), color = "gray80", linetype = "dashed") +
  geom_col(position = "stack") +
  scale_fill_manual(
    values = c(custom_class_colors, Majority = "#999999"),
    name = NULL,
    labels = ec_labels
  ) +
  scale_y_continuous(breaks = seq(0, 25, by = 5), limits = c(0, 25)) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = "right"
  )

# === Combine plots A + B ===
final_plot <- annotation_plot + heatmap + plot_layout(ncol=2, widths=c(0.5, 4))

ggarranged_combined_plot <- ggarrange(
  final_plot + theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10)),                 
  top1_barplot_with_majority + theme(plot.margin = margin(t = 30, r = 0, b = 10, l = 50)),        
  labels = c("A", "B"),              
  label.x = 0.01,                    
  label.y = 1,                     
  ncol = 1, nrow = 2,                
  heights = c(1.3, 1),               
  common.legend = FALSE,
  font.label = list(size = 16, face = "bold")  
)

print(ggarranged_combined_plot)

# === Save final figure ===
ggsave(
  filename = "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/CaseStudy/casestudyplot3.png",    
  plot = ggarranged_combined_plot,            
  width = 10,                      
  height = 15,                      
  dpi = 300,   
  bg= "white"
)

# === Save Majority Vote Table (only if included) ===
if (include_majority_vote) {
  majority_table <- majority_preds %>%
    select(drug, true_EC = drug_EC, 
           majority_pred_top1 = pred_top1, 
           majority_pred_top5 = pred_top5, 
           top1_hit_type, 
           top5_hit_type)
  
  write_csv(
    majority_table, 
    "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/CaseStudy/majority_vote_results3.csv"
  )
}
