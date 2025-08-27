library(tidyverse)
library(scales)
library(cowplot)
library(paletteer)

## ----- 1) NEEDLE DATA (pick the samples you want) -----
kept_indv <- c(
  "CHES_P12002_145", "CHES_P12002_146", "CHES_P12002_147", "CHES_P12002_148", 
  "CHES_P12002_150", "CHES_P12002_152", "CHES_P12002_153", "CHFI_P12002_122", 
  "CHFI_P12002_123", "CHFI_P12002_125", "CHFI_P12002_126", "CHFI_P12002_128", 
  "CHFI_P12002_129", "CHFI_P12002_130", "CHFI_P12002_131", "CHFI_P12002_132", 
  "CHSK_P12002_101", "CHSK_P12002_102", "CHSK_P12002_103", "CHSK_P12002_104", 
  "CHSK_P12002_105", "CHSK_P12002_106", "CHSK_P12002_107", "CHSK_P12002_108", 
  "CHSK_P12002_109", "CHSK_P12002_110", "CHSK_P14052_101", "CHSK_P14052_102", 
  "CHST_P12002_112", "CHST_P12002_113", "CHST_P12002_114", "CHST_P12002_115", 
  "CHST_P12002_116", "CHST_P12002_121", "COES_P12002_175", "COES_P12002_176", 
  "COES_P12002_177", "COES_P12002_179", "COES_P12002_182", "COES_P12002_183", 
  "COES_P12002_184", "COGE_P12002_133", "COGE_P12002_135", "COGE_P12002_136", 
  "COGE_P12002_137", "COGE_P12002_138", "COGE_P12002_139", "COLI_P12002_164", 
  "COLI_P12002_165", "COLI_P12002_166", "COLI_P12002_167", "COLI_P12002_168", 
  "COLI_P12002_169", "COLI_P12002_170", "COLI_P12002_172", "COLI_P12002_173", 
  "COLI_P12002_174", "COSK_P12002_154", "COSK_P12002_155", "COSK_P12002_156", 
  "COSK_P12002_157", "COSK_P12002_158", "COSK_P12002_159", "COSK_P12002_161", 
  "COSK_P12002_162", "COSK_P12002_163", "COSK_P14052_109"
)

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

rank <- "genus"

bamqc <- read_delim("data/stammerula2025/stats/qualimap-summary-stats.csv", delim = ",")

df_imiss2 <- read_delim("data/stammerula2025/stats/stam-stam.filt.bial.snps.imiss", delim = "\t")

imiss_df <- df_imiss2 %>%
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

imiss_for_plot <- df %>% 
  filter(!REGION == "scot")

# Shared order (highest missingness first)
ord <- imiss_for_plot %>%
  arrange(desc(F_MISS)) %>%
  pull(sample)

imiss_for_plot <- imiss_for_plot %>%
  mutate(sample = factor(sample, levels = ord))

## ----- 2) STACKED BAR DATA (species) -----

df <- read_delim(paste0("data/KRAKEN/rel_abund_0_2/ra_", rank,".csv"),
                 delim = ",", col_names = T)


df_plot <- df %>% 
  filter(taxon != "Tephritis") %>% 
  mutate(rank = if_else(rel_abund_perc < 0.01, "Other", taxon))


# Now collapse all the "Other" rows per sample
df_plot2 <- df_plot %>%
  filter(!(Sample_id %in% c("pt_042_cell1", "pt_042_cell2", "pt_042_cell3"))) %>% 
  mutate(Sample_id = factor(Sample_id, levels = ord))

## ----- 3) MAKE PLOTS (vertical orientation; no coord_flip) -----
# Top: needle plot
p_needles <- ggplot(imiss_for_plot, aes(x = sample, y = F_MISS)) +
  geom_segment(aes(xend = sample, y = 0, yend = F_MISS), linewidth = 1, color = "grey55") +
  geom_point(size = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(x = NULL, y = "Missingness (F_MISS)") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_numb <- ggplot(imiss_for_plot, aes(x = sample, y = n_mapped_reads)) +
  geom_segment(aes(xend = sample, y = 0, yend = F_MISS), linewidth = 1, color = "grey55") +
  geom_point(size = 1) +
  #scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(x = NULL, y = "Missingness (F_MISS)") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# compute scaling factor between F_MISS and log10_aligned
f_max <- max(imiss_for_plot$F_MISS, na.rm = TRUE)
a_max <- max(imiss_for_plot$n_reads, na.rm = TRUE)
ratio <- f_max / a_max

p_top <- ggplot(imiss_for_plot, aes(x = sample)) +
  # F_MISS line (left axis, unscaled)
  geom_line(aes(y = F_MISS, color = "F_MISS"), linewidth = 1) +
  geom_point(aes(y = F_MISS, color = "F_MISS"), size = 0.8) +
  # log10_aligned line (scaled to match F_MISS axis for plotting)
  geom_line(aes(y = n_reads * ratio, color = "number of reads"), linewidth = 0.6) +
  geom_point(aes(y = n_reads * ratio, color = "numberreads"), size = 1.5) +
  # dual y-axes
  scale_y_continuous(
    name = "Missingness (F_MISS)",
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    sec.axis = sec_axis(~ . / ratio, name = "number of decon. reads",
                        labels = label_number(scale_cut = cut_short_scale()))
  ) +
  labs(x = NULL, color = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.line.y.left = element_line(color = "black"),
    axis.line.y.right = element_line(color = "black"),
    axis.ticks.y.left = element_line(color = "black"),
    axis.ticks.y.right = element_line(color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top"
  )

colors <- paletteer_d("colorBlindness::SteppedSequential5Steps")

# Bottom: stacked relative-abundance bars
p_bars <- ggplot(df_plot2, aes(x = Sample_id, y = rel_abund_perc, fill = rank)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL, y = "Relative abundance (bacteria)", fill = rank) +
  theme_minimal(base_size = 12) +
  #scale_fill_brewer(palette = "Paired") +
  #scale_fill_manual(values = colors) +
  scale_fill_paletteer_d("ggthemes::Classic_20",
                         guide = guide_legend(nrow = 3, byrow = TRUE)) +
  #scale_fill_paletteer_d("ggthemes::Tableau_20") +
  #scale_fill_paletteer_d("ggsci::category20b_d3") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 8),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "inches")
  )

## ----- 4) COMBINE WITH COWPLOT -----
combo <- plot_grid(
  #p_needles + theme(legend.position = "none"),
  #p_numb + theme(legend.position = "none"),
  p_top + theme(legend.position = "none"),
  p_bars,
  ncol = 1, align = "v", axis = "lr",
  rel_heights = c(0.8, 2.2)   # adjust to taste
)

combo

# Save
ggsave(paste0("00_thesis-fig/k2-fixed/trio_rel_abund_all",rank,".png"), combo, width = 8, height = 6, dpi = 300, units = "in")
ggsave(paste0("00_thesis-fig/k2-fixed/trio_rel_abund_all",rank,".pdf"), combo, width = 8, height = 6, device = cairo_pdf)
ggsave(paste0("00_thesis-fig/k2-fixed/trio_rel_abund_all",rank,".svg"), combo, width = 12, height = 8, device = "svg", units = "in", dpi = 300)


ggsave(paste0("00_thesis-fig/k2-fixed/trio_rel_abund_all",rank,"_top.svg"), p_top, width = 12, height = 8, dpi = 300, units = "in")
ggsave(paste0("00_thesis-fig/k2-fixed/trio_rel_abund_all",rank,"_bars.svg"), p_bars, width = 12, height = 8, dpi = 300, units = "in")
