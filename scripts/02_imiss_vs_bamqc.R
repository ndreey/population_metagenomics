library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(GGally)

# ----- PARAMETERS -----
vcf <- "stam.filt.bial.snps"
output_dir <- paste0("00_thesis-fig/fig1/",vcf)
imiss_filt <- 0.7
device_type <- "svg"  # "svg" or "png"
plot_width <- 10
plot_height <- 5
plot_dpi <- 300
plot_units <- "in"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----- FUNCTION TO SAVE PLOTS -----

save_plot <- function(plot_obj, filename, width = plot_width, height = plot_height,
                      dpi = plot_dpi, device = device_type, units = plot_units, bg = "white") {
  full_path <- file.path(output_dir, paste0(filename, ".", device))
  ggsave(filename = full_path, plot = plot_obj, width = width, height = height,
         dpi = dpi, device = device, units = units, bg = bg)
}



# ----- Load in Data -----
df_imiss_1 <- read_delim("data/stammerula2025/stats/stam-stam.bial.snps.imiss", delim = "\t")
df_imiss_2 <- read_delim("data/stammerula2025/stats/stam-stam.filt.bial.snps.imiss", delim = "\t")
bamqc <- read_delim("data/stammerula2025/stats/qualimap-summary-stats.csv", delim = ",")

df_imiss <- df_imiss_2

# ----- METADATA -----
metadata <- tribble(
  ~POP,  ~HP, ~REGION, ~regHP,
  "CHES", "CH", "east", "eastCH",
  "CHFI", "CH", "east", "eastCH",
  "CHSK", "CH", "west", "westCH",
  "CHST", "CH", "west", "westCH",
  "COES", "CO", "east", "eastCO",
  "COGE", "CO", "west", "westCO",
  "COLI", "CO", "east", "eastCO",
  "COSK", "CO", "west", "westCO",
  "CHSC", "SC", "scot", "scotHP",
  "CPSC", "SC", "scot", "scotHP",
)

# ----- COLOR SCHEMES -----

dark2 <- brewer.pal(8, "Dark2")

p_colors <- list(
  HP = c(
    "CH" = "orchid",
    "CO" = "palegreen3",
    "SC" = "slateblue"),
  REGION = c(
    "west" = "wheat",
    "east" = "wheat4",
    "scot" = "steelblue1"),
  regHP = c(
    "westCH" = "orchid1", 
    "eastCH" = "orchid4", 
    "westCO" = "palegreen1", 
    "eastCO" = "palegreen4",
    "scotHP" = "skyblue3"),
  POP = c(
    "CHES" = dark2[1],  
    "CHFI" = dark2[2], 
    "CHSC" = "cadetblue", 
    "CHSK" = dark2[3],
    "CHST" = dark2[4],  
    "COES" = dark2[5],  
    "COGE" = dark2[6],  
    "COLI" = dark2[7],  
    "COSK" = dark2[8], 
    "CPSC" = "cadetblue1")
)
# ----- PROCESSING -----

imiss_df <- df_imiss %>%
  mutate(
    POP = str_extract(INDV, "^[^_]+")
  ) %>%
  left_join(metadata, by = "POP") %>%
  mutate(
    REGION = factor(REGION, levels = c("west", "east", "scot")),
    HP = factor(HP, levels = c("CH", "CO", "SC")),
    regHP = factor(regHP, levels = c("westCH", "eastCH", "westCO", "eastCO", "scotHP"))
  )

bamqc_df <- bamqc %>%
  mutate(
    POP = str_extract(sample, "^[^_]+")
  ) %>%
  left_join(metadata, by = "POP") %>%
  mutate(
    REGION = factor(REGION, levels = c("west", "east", "scot")),
    HP = factor(HP, levels = c("CH", "CO", "SC")),
    regHP = factor(regHP, levels = c("westCH", "eastCH", "westCO", "eastCO", "scotHP"))
  )

# Add imiss_df$F_MISS to bamqc_df
df <- bamqc_df %>%
  left_join(imiss_df %>% select(INDV, F_MISS), by = c("sample" = "INDV"))

df <- df %>%
  mutate(
    log10_reads = log10(n_reads),
    log10_aligned = log10(n_mapped_reads),
  )

df <- df %>% 
  filter(!REGION == "scot")

# ---- PLOTS -----

# F_MISS vs aligned reads
plot_aligned <- df %>%
  ggplot(aes(x = log10_aligned, y = F_MISS)) +
  geom_point(color = "black", size = 1.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Aligned Reads vs Missingness",
    x = "log10(Aligned Reads)",
    y = "Individual Missingness (F_MISS)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )


# F_MISS vs Total Reads
plot_tot_reads <- df %>% 
  ggplot(aes(x = n_reads, y = F_MISS)) +
  geom_point(color = "black", size = 1.5) +
  #geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = pretty_breaks(n=10), 
                     labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "Total Reads vs Missingness",
    x = "Total Reads",
    y = "Individual Missingness (F_MISS)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )

# Mean Coverage vs F_MISS
plot_mean_cov <- df %>%
  ggplot(aes(x = mean.Cov, y = F_MISS)) +
  geom_point(color = "black", size = 1.5) +
  #geom_smooth(method = "lm", color = "blue", se = FALSE, linewidth = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = pretty_breaks(n=20), limit = c(0,10)) +
  labs(
    title = "Mean Coverage vs Missingness",
    x = "Mean Coverage",
    y = "Individual Missingness (F_MISS)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank()
  )


# ===== imiss GG PLOTS =====

# Compare with GGPAIRs
p_pair <- df %>%
  filter(!REGION == "scot") %>%
  select(F_MISS, log10_aligned, mean.Cov, GC, log10_reads, HP) %>%
  ggpairs(
    columns = 1:5,
    mapping = aes(color = HP),  #
    lower = list(
      continuous = wrap("points", alpha = 0.4, size = 1, color = "black")
    ),
    upper = list(
      continuous = wrap("cor", size = 2.5)
    ),
    diag = list(
      continuous = wrap("densityDiag", alpha = 0.6)
    ),
    legend = grab_legend(
      ggplot(df, aes(x = F_MISS, fill = HP)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = p_colors$HP)
    )
  ) +
  scale_fill_manual(values = p_colors$HP) +  # Apply to the full plot
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )





# ----- Flexible Parameters -----
x_var <- "log10_aligned"  # Choose: "mean.Cov", "log10_aligned", "GC", etc.
x_label <- "log10(Aligned Reads)"
x_threshold <- 3.8
f_miss_threshold <- 0.7

# ----- Identify Kept Individuals -----
kept_indv <- df %>%
  filter(!REGION == "scot") %>%
  filter(F_MISS <= f_miss_threshold, .data[[x_var]] >= x_threshold) %>%
  pull(sample)

# ----- Summary per POP -----
pop_summary_all <- df %>%
  filter(!REGION == "scot") %>%
  group_by(POP) %>%
  summarise(
    n_total = n(),
    n_kept = sum(F_MISS <= f_miss_threshold & .data[[x_var]] >= x_threshold),
    .groups = "drop"
  )

pop_lines <- pop_summary_all %>%
  mutate(line = paste0(POP, ": ", n_kept, " (", round(100 * n_kept / n_total, 1), "%)")) %>%
  pull(line) %>%
  paste(collapse = "\n")

total_n <- sum(pop_summary_all$n_total)
total_kept <- sum(pop_summary_all$n_kept)
total_frac <- round(100 * total_kept / total_n, 1)

final_label <- paste0("Total kept: ", total_kept, " / ", total_n, " (", total_frac, "%)")

# ----- Plot -----
imiss_xvar_all_p <- df %>%
  filter(!REGION == "scot") %>%
  ggplot(aes(x = .data[[x_var]], y = F_MISS)) +
  geom_point(size = 1, alpha = 0.6, color = "grey40") +
  geom_smooth(method = "gam",
              formula = y ~ s(x, bs = "cs"),
              se = FALSE,
              color = "black",
              linetype = "solid",
              linewidth = 0.5) +
  geom_hline(yintercept = f_miss_threshold, linetype = "dashed",
             color = "firebrick4", alpha = 0.6, linewidth = 0.3) +
  geom_vline(xintercept = x_threshold, linetype = "dashed",
             color = "firebrick2", alpha = 0.6, linewidth = 0.3) +
  geom_text(x = x_threshold, y = 0.9, label = "X threshold",
            color = "black", size = 2, vjust = -1, angle = 90) +
  geom_text(x = mean(df[[x_var]], na.rm = TRUE), y = f_miss_threshold,
            label = "F_MISS threshold",
            color = "black", size = 2, vjust = -1) +
  annotate("text", x = mean(df[[x_var]], na.rm = TRUE), y = 0.95, label = final_label,
           size = 5, color = "black", hjust = 0) +
  annotate("label", x = quantile(df[[x_var]], 0.95), y = 0.5, label = pop_lines,
           size = 3, color = "black", hjust = 1) +
  geom_text_repel(data = df %>% filter(sample %in% kept_indv),
                  aes(label = sample),
                  size = 1,
                  max.overlaps = 20,
                  segment.size = 0.2,
                  min.segment.length = 0,
                  segment.color = "grey60",
                  box.padding = 0.2,
                  point.padding = 0.1,
                  show.legend = FALSE) +
  scale_y_continuous(breaks = pretty_breaks(n = 20),
                     labels = label_number(scale_cut = cut_short_scale())) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 20),
    labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),
    legend.position = "right"
  ) +
  labs(
    x = x_label,
    y = "Fraction Missing (F_MISS)",
    title = paste(x_label, "vs. Missingness"),
    subtitle = paste0(x_var, " threshold: ", x_threshold,
                      "\nF_MISS threshold: ", f_miss_threshold)
  )


# ----- Create Subset of Passing Individuals -----
df_passed <- df %>%
  filter(F_MISS <= f_miss_threshold, .data[[x_var]] >= x_threshold)

# ----- Pair Plot for Filtered Samples -----
p_pair_passed <- df_passed %>%
  filter(!REGION == "scot") %>%
  select(F_MISS, log10_aligned, mean.Cov, GC, log10_reads, HP) %>%
  ggpairs(
    columns = 1:5,
    mapping = aes(color = HP),  #
    lower = list(
      continuous = wrap("points", alpha = 0.4, size = 1, color = "black")
    ),
    upper = list(
      continuous = wrap("cor", size = 2.5)
    ),
    diag = list(
      continuous = wrap("densityDiag", alpha = 0.6)
    ),
    legend = grab_legend(
      ggplot(df_passed, aes(x = F_MISS, fill = HP)) +
        geom_density(alpha = 0.5) +
        scale_fill_manual(values = p_colors$HP)
    )
  ) +
  scale_fill_manual(values = p_colors$HP) + 
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


# ----- Save the plots -----
save_plot(plot_aligned, filename = "log10_aligned_vs_F_MISS")
save_plot(plot_tot_reads, filename = "totReads_vs_F_MISS")
save_plot(plot_mean_cov, filename = "cov_vs_F_MISS")
save_plot(p_pair, filename = "pairplot_HP_vs_qualimap")
save_plot(imiss_xvar_all_p, paste0("imiss_vs_", x_var, "_summary"))
save_plot(p_pair_passed, paste0("pairplot_HP_vs_qualimap_passed_", x_var))


# Save list of samples to .txt file
write_lines(df_passed$sample, "data/stammerula2025/stats/samples_passed_qualimap.txt")





# Plot F_MISS vs mean.Cov
df_passed %>%
  ggplot(aes(x = mean.Cov, y = F_MISS)) +
  geom_point(size = 1, alpha = 0.6, color = "grey40") +
  xlim(0, 10) +
  scale_x
  









