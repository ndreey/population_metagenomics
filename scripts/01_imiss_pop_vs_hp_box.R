library(tidyverse)
library(RColorBrewer)
library(scales)

# ----- PARAMETERS -----
vcf <- "stam.filt.bial.snps"
output_dir <- paste0("00_thesis-fig/fig1/",vcf)
imiss_filt <- 0.7
device_type <- "svg"  # "svg" or "png"
plot_width <- 7
plot_height <- 6
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

df_imiss <- df_imiss_2

# ----- METADATA -----
# Metadata
metadata <- tribble(
  ~POP,  ~HP, ~REGION, ~regHP, ~contact, ~color, ~shape,
  "CHES", "CH", "east", "eastCH", "sympatric", "#BCBDDC", 23,
  "CHFI", "CH", "east", "eastCH", "allopatric", "#6A51A3", 23,
  "CHSK", "CH", "west", "westCH", "sympatric", "#9E9AC8", 21,
  "CHST", "CH", "west", "westCH", "allopatric", "#4A1486", 21,
  "COES", "CO", "east", "eastCO", "sympatric", "#C7E9C0", 23,
  "COGE", "CO", "west", "westCO", "allopatric", "#00441B", 21,
  "COLI", "CO", "east", "eastCO", "allopatric", "#238B45", 23,
  "COSK", "CO", "west", "westCO", "sympatric", "#74C476", 21,
  "CHSC", "CH", "scot", "scotHP", "sympatric", "cadetblue", 20,
  "CPSC", "CP", "scot", "scotHP", "sympatric", "cadetblue1", 20
)

p_colors <- list(
  HP = c(
    "CH" = "#4A1486",
    "CO" = "#00441B"),
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
    "CHES" = "#BCBDDC",  
    "CHFI" = "#6A51A3", 
    "CHSC" = "cadetblue", 
    "CHSK" = "#9E9AC8",
    "CHST" = "#4A1486",  
    "COES" = "#C7E9C0",  
    "COGE" = "#00441B",  
    "COLI" = "#238B45",  
    "COSK" = "#74C476", 
    "CPSC" = "cadetblue1")
)

p_shape <- list(
  POP = c(
    "CHES" = 23,  
    "CHFI" = 23, 
    "CHSC" = 20, 
    "CHSK" = 21,
    "CHST" = 21,  
    "COES" = 23,  
    "COGE" = 21,  
    "COLI" = 23,  
    "COSK" = 21, 
    "CPSC" = 20)
)

# ----- PROCESSING -----

imiss_df <- df_imiss %>%
  mutate(
    POP = str_extract(INDV, "^[^_]+")
  ) %>%
  left_join(metadata, by = "POP") %>%
  mutate(
    REGION = factor(REGION, levels = c("west", "east", "scot")),
    HP = factor(HP, levels = c("CH", "CO")),
    regHP = factor(regHP, levels = c("westCH", "eastCH", "westCO", "eastCO", "scotHP"))
  )

# ----- PLOTTING -----


# Boxplot of F_MISS by POP, faceted by HP (excluding scots)
imiss_p_box_POPxHP <- imiss_df %>%
  filter(REGION != "scot") %>%
  ggplot(aes(x = POP, y = F_MISS, fill = POP)) +
  geom_boxplot(alpha = 0.8, outliers = F) +
  geom_jitter(alpha = 0.4, width = 0.1, shape = 21, fill = "black", color = "black", size = 1.5) +
  scale_fill_manual(values = p_colors$POP) +
  scale_y_continuous(breaks = pretty_breaks(n=15)) +
  facet_wrap(~HP, scales = "free_x") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.grid   = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2),  # bold border
    strip.text = element_text(size = 16),
    strip.background = element_blank()
  ) +
  labs(
    x = "Population",
    y = "Individual Missingness (F_MISS)"
  )


save_plot(imiss_p_box_POPxHP, "UPDATED_imiss-box_POPxHP")


bamqc <- read_delim("data/stammerula2025/stats/qualimap-summary-stats.csv", delim = ",")


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

df

x_threshold <- 3.8
f_miss_threshold <- 0.7

# Plot log10_aligned vs F_MISS
imis_vs_aligned <- df %>% 
  ggplot(aes(x=log10_aligned, y=F_MISS)) +
  geom_point(size=3, alpha=0.7) +
  #scale_fill_manual(values = p_colors$POP) +
  #scale_shape_manual(values = p_shape$POP) +
  geom_hline(yintercept = f_miss_threshold, linetype = "dashed",
             color = "black", alpha = 0.6, linewidth = 0.3) +
  geom_vline(xintercept = x_threshold, linetype = "dashed",
             color = "black", alpha = 0.6, linewidth = 0.3) +
  scale_y_continuous(breaks = pretty_breaks(n = 15),
                     labels = label_number(scale_cut = cut_short_scale())) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 15),
    labels = label_number(scale_cut = cut_short_scale())) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    panel.grid   = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1.2),  # bold border
    strip.text = element_text(size = 16),
    strip.background = element_blank()
  )
# Save as svg
ggsave("00_thesis-fig/ALIGNED_VS_FMISS.svg", 
       imis_vs_aligned, width = 10, height = 6)
