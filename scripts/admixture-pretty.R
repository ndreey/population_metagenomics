library(tidyverse)
library(scales)
library(cowplot)
library(paletteer)

# --- Load metadata ---

metadata <- tribble(
  ~POP,  ~HP, ~REGION, ~regHP, ~contact,
  "CHES", "CH", "east", "eastCH", "sympatric",
  "CHFI", "CH", "east", "eastCH", "allopatric",
  "CHSK", "CH", "west", "westCH", "sympatric",
  "CHST", "CH", "west", "westCH", "allopatric",
  "COES", "CO", "east", "eastCO", "sympatric",
  "COGE", "CO", "west", "westCO", "allopatric",
  "COLI", "CO", "east", "eastCO", "allopatric",
  "COSK", "CO", "west", "westCO", "sympatric",
  "CHSC", "CH", "scot", "scotHP", "sympatric",
  "CPSC", "CP", "scot", "scotHP", "sympatric"
)

POP_ORDER <- c("CHSK", "CHST", "CHFI", "CHES", "COES", "COLI", "COGE", "COSK")

# Load color palettes
green <- rev(paletteer_d("ggsci::light_green_material"))[-1]
purple <- rev(paletteer_d("ggsci::purple_material"))[-1]

# Define color schemes by K
color_schemes <- list(
  `2` = c("Q1" = "mediumpurple4", "Q2" = "mediumpurple3"),
  `3` = c("Q1" = "mediumpurple4", "Q2" = "mediumpurple1", "Q3" = "palegreen3"),
  `6` = c("Q1" = "palegreen3", "Q2" = "mediumpurple3", "Q3" = "mediumpurple4", 
        "Q4" = "palegreen4", "Q5" = "cadetblue", "Q6" = "cadetblue1")
)

# Loop over K values
K_values <- c(2, 3, 6)
all_plots <- list()

for (K in K_values) {
  
  # Load Q matrix
  df_q <- read_delim(
    paste0("data/stammerula2025/admixture/stam.passed.no-sc.no-outlier.", K, ".Q"),
    delim = " ", col_names = FALSE
  )
  colnames(df_q) <- paste0("Q", 1:K)
  
  # Load sample IDs
  meta_df <- read_delim("data/stammerula2025/admix.fam.samples.txt",
                        delim = "\t", col_names = "INDV") %>%
    mutate(POP = str_extract(INDV, "^[^_]+"))
  
  # Combine Q with metadata
  df_q$INDV <- meta_df$INDV
  long_df <- df_q %>%
    pivot_longer(cols = starts_with("Q"), names_to = "Cluster", values_to = "Proportion") %>%
    left_join(meta_df, by = "INDV") %>%
    left_join(metadata, by = "POP") %>%
    group_by(POP) %>%
    arrange(POP, INDV) %>%
    ungroup()
  
  long_df$POP <- factor(long_df$POP, levels = POP_ORDER)
  
  ordered_ids <- long_df %>%
    filter(Cluster == "Q2") %>%
    arrange(factor(POP, levels = POP_ORDER), desc(Proportion)) %>%
    pull(INDV) %>%
    unique()
  
  # 2. Update INDV to be a factor in that order
  long_df <- long_df %>%
    mutate(INDV = factor(INDV, levels = ordered_ids))
  
  # Plot
  p_main <- ggplot(long_df, aes(x = INDV, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.99) +
    facet_grid(~POP, scales = "free_x", space = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = color_schemes[[as.character(K)]]) +
    #scale_fill_brewer(palette = "Paired", direction = -1) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(size = 14),
      legend.position = "none"
    )
  
  # Right-side K label
  label_plot <- ggdraw() +
    draw_label(paste0("K = ", K), angle = 270, size = 12)
  
  # Combine plot and label
  combined <- plot_grid(p_main, label_plot, ncol = 2, rel_widths = c(1, 0.03))
  all_plots[[as.character(K)]] <- combined
}

# Stack all plots vertically
final_plot <- plot_grid(plotlist = all_plots, ncol = 1)
final_plot
# Save
ggsave("00_thesis-fig/admixture/admixture_K2-3_panel-pretty.svg",
       final_plot, width = 14, height = 2.2 * length(K_values), dpi = 300)



## Plot cv error plot
cv_error <- read_delim("data/stammerula2025/admixture/cv_error.tsv", 
                       delim = "\t", col_names = c("K", "cv_error"))

cv_error <- cv_error %>%
  arrange(K) %>% 
  filter(K <= 11)

cv_plot <- ggplot(cv_error, aes(x = K, y = cv_error)) +
  geom_line(color = "blue") +
  geom_point(size = 2, color = "blue") +
  labs(x = "K", y = "Cross-validation error") +
  scale_x_continuous(breaks = pretty_breaks(n = 20)) +
  scale_y_continuous(breaks = pretty_breaks(n=10)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave("00_thesis-fig/admixture/admixture_cv_error.png",
       cv_plot, width = 5, height = 4, dpi = 300)

