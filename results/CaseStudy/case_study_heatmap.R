# === Parameters ===
include_majority_vote <- TRUE   # set to FALSE to exclude majority vote from plots

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

# Set seed for reproducible ordering of tied predictions
set.seed(42)

# READ DATAFRAME HERE
df <- read_csv("/Users/josefinaarcagni/Downloads/merged_output_case3.csv")

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

# === Hierarchical tie-breaking for majority vote ===
select_winner_hierarchical <- function(rank_votes) {
  # rank_votes is a list where each EC has a vector of counts [rank1_count, rank2_count, ...]
  
  if (length(rank_votes) == 0) return(NA)
  
  # Convert to data frame for easier manipulation
  vote_df <- data.frame(
    ec = names(rank_votes),
    stringsAsFactors = FALSE
  )
  
  # Add rank counts
  max_ranks <- max(sapply(rank_votes, length))
  for (i in 1:max_ranks) {
    vote_df[[paste0("rank", i)]] <- sapply(rank_votes, function(x) {
      if (length(x) >= i) x[i] else 0
    })
  }
  
  # Sort by rank1 (descending), then rank2, then rank3, etc.
  rank_cols <- paste0("rank", 1:max_ranks)
  vote_df <- vote_df %>%
    arrange(across(all_of(rank_cols), desc))
  
  # Get the top scorer(s)
  top_row <- vote_df[1, ]
  
  # Check if there are ties (same counts across all ranks)
  if (nrow(vote_df) > 1) {
    is_tie <- TRUE
    for (col in rank_cols) {
      if (top_row[[col]] != vote_df[2, col]) {
        is_tie <- FALSE
        break
      }
    }
    
    if (is_tie) {
      # Find all ECs with identical rank counts
      tied_ecs <- c(top_row$ec)
      for (i in 2:nrow(vote_df)) {
        same <- TRUE
        for (col in rank_cols) {
          if (vote_df[i, col] != top_row[[col]]) {
            same <- FALSE
            break
          }
        }
        if (same) {
          tied_ecs <- c(tied_ecs, vote_df$ec[i])
        } else {
          break
        }
      }
      # Sort and concatenate
      tied_ecs <- sort(tied_ecs)
      return(paste(tied_ecs, collapse = "|"))
    }
  }
  
  return(top_row$ec)
}

# === Helper function to check if prediction matches true EC and return fractional hit ===
check_hit_with_pipes_fractional <- function(prediction, true_prefix) {
  if (is.na(prediction) || prediction == "") return(0)
  
  # Split by pipe to handle tied predictions (same rank)
  pred_options <- str_split(prediction, "\\|")[[1]]
  pred_options <- str_trim(pred_options)
  
  # Count how many match
  num_matches <- sum(pred_options == true_prefix)
  
  # Return fractional hit: 1/n where n is total number of tied predictions
  if (num_matches > 0) {
    return(1 / length(pred_options))
  } else {
    return(0)
  }
}

# === Prepare data ===
ec_methods <- c("E-zyme1","E-zyme2", "BridgIT", "SelenzymeRF","SIMMER", "Theia" ,"BEC-Pred","CLAIRE")
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
    class = paste0(class),
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

# === Majority vote with Top-1 and Top-5 (with fractional hits) ===
if (include_majority_vote) {
  majority_methods <- c("Theia", "BEC-Pred", "SIMMER", "SelenzymeRF")
  
  majority_preds <- df %>%
    select(drug, drug_EC, all_of(majority_methods)) %>%
    rowwise() %>%
    mutate(
      true_prefix = paste(str_split(drug_EC, "\\.")[[1]][1:3], collapse="."),
      
      # Top-1 Majority vote: only use rank 1
      pred_majority_top1 = {
        rank_votes <- list()
        
        for (m in majority_methods) {
          pred_str <- get(m)
          if (!is.na(pred_str) && pred_str != "") {
            # Get only the first rank (Top-1)
            preds <- str_split(pred_str, ";")[[1]]
            if (length(preds) >= 1) {
              pipes <- str_split(preds[[1]], "\\|")[[1]]
              pipes <- str_trim(pipes)
              pipes <- pipes[pipes != ""]
              
              for (p in pipes) {
                subsub <- paste(str_split(p, "\\.")[[1]][1:3], collapse=".")
                if (is.null(rank_votes[[subsub]])) {
                  rank_votes[[subsub]] <- 1
                } else {
                  rank_votes[[subsub]] <- rank_votes[[subsub]] + 1
                }
              }
            }
          }
        }
        
        # Select winner based on votes
        if (length(rank_votes) == 0) {
          NA
        } else {
          max_votes <- max(unlist(rank_votes))
          winners <- names(rank_votes)[unlist(rank_votes) == max_votes]
          if (length(winners) > 1) {
            paste(sort(winners), collapse = "|")
          } else {
            winners[1]
          }
        }
      },
      
      # Top-5 Majority vote: use ranks 1-5 with hierarchical tie-breaking
      pred_majority_top5 = {
        rank_votes <- list()
        
        for (m in majority_methods) {
          pred_str <- get(m)
          if (!is.na(pred_str) && pred_str != "") {
            preds <- str_split(pred_str, ";")[[1]]
            preds <- preds[1:min(5, length(preds))]
            
            for (rank_idx in seq_along(preds)) {
              pipes <- str_split(preds[[rank_idx]], "\\|")[[1]]
              pipes <- str_trim(pipes)
              pipes <- pipes[pipes != ""]
              
              for (p in pipes) {
                subsub <- paste(str_split(p, "\\.")[[1]][1:3], collapse=".")
                
                if (is.null(rank_votes[[subsub]])) {
                  rank_votes[[subsub]] <- rep(0, 5)
                }
                
                rank_votes[[subsub]][rank_idx] <- rank_votes[[subsub]][rank_idx] + 1
              }
            }
          }
        }
        
        select_winner_hierarchical(rank_votes)
      },
      
      # Check matches with fractional hits
      majority_fractional_hit_top1 = check_hit_with_pipes_fractional(pred_majority_top1, true_prefix),
      majority_fractional_hit_top5 = check_hit_with_pipes_fractional(pred_majority_top5, true_prefix),
      
      # Classify hit types (for display)
      majority_hit_type_top1 = case_when(
        is.na(pred_majority_top1) | pred_majority_top1 == "" ~ "No prediction",
        majority_fractional_hit_top1 > 0 ~ "Top 1",
        TRUE ~ "No hit"
      ),
      
      majority_hit_type_top5 = case_when(
        is.na(pred_majority_top5) | pred_majority_top5 == "" ~ "No prediction",
        majority_fractional_hit_top5 > 0 ~ "Top 5",
        TRUE ~ "No hit"
      )
    ) %>% ungroup()
  
  # Individual method Top-1 hits
  top1_hits_only <- heatmap_data %>%
    filter(hit_type == "Top 1") %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    count(method, class) %>%
    mutate(HitType = "Top 1")
  
  # Majority vote Top-1 hits (with fractional counting)
  majority_top1 <- majority_preds %>%
    filter(majority_fractional_hit_top1 > 0) %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    group_by(method = "Top-1", class) %>%
    summarise(n = sum(majority_fractional_hit_top1), .groups = "drop") %>%
    mutate(HitType = "Majority Top-1")
  
  # Majority vote Top-5 hits (with fractional counting)
  majority_top5 <- majority_preds %>%
    filter(majority_fractional_hit_top5 > 0) %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    group_by(method = "Top-5", class) %>%
    summarise(n = sum(majority_fractional_hit_top5), .groups = "drop") %>%
    mutate(HitType = "Majority Top-5")
  
  top1_plus_majority <- bind_rows(top1_hits_only, majority_top1, majority_top5)
  
} else {
  top1_plus_majority <- heatmap_data %>%
    filter(hit_type == "Top 1") %>%
    left_join(tree_data_ordered %>% select(drug, class), by="drug") %>%
    count(method, class) %>%
    mutate(HitType = "Top 1")
}

# Set factor levels depending on inclusion
if (include_majority_vote) {
  top1_plus_majority <- top1_plus_majority %>%
    mutate(method = factor(method, levels = c(ec_methods, "Top-1", "Top-5")))
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
  scale_y_continuous(breaks = seq(0, 27, by = 5), limits = c(0, 27)) +
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
  filename = "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/CaseStudy/casestudyplot_final.jpg",    
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
           majority_prediction_top1 = pred_majority_top1,
           majority_fractional_hit_top1,
           majority_hit_type_top1,
           majority_prediction_top5 = pred_majority_top5,
           majority_fractional_hit_top5, 
           majority_hit_type_top5)
  
  write_csv(
    majority_table, 
    "/Users/josefinaarcagni/Documents/ECMethods/FinalGraphs/CaseStudy/majority_vote_results_final.csv"
  )
}