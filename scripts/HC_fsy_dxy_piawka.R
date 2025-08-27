library(tidyverse)
library(scales)
library(paletteer)
library(shadowtext)
library(cowplot)
library(pheatmap)

# ----- PARAMETERS -----
input_tsv <- "data/stammerula2025/fst/POP-passed.all-sites.no-sc.hudson.tsv"
output_dir <- "00_thesis-fig/piawka"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ----- LOAD FULL DATA -----
df <- read_tsv(input_tsv, col_types = cols())

# ---------------------------------------------------
fst_df <- df %>%
  filter(metric == "Fst_HUD" & !is.na(pop2) & pop2 != ".") %>%
  select(pop1, pop2, value) %>%
  mutate(value = as.numeric(value))

# Make symmetric Fst matrix with 0 diagonal
fst_mat <- fst_df %>%
  # Add the reversed pairs to ensure symmetry
  bind_rows(fst_df %>% rename(pop1 = pop2, pop2 = pop1)) %>%
  # Add zeros for diagonal
  bind_rows(tibble(pop1 = unique(c(fst_df$pop1, fst_df$pop2)),
                   pop2 = unique(c(fst_df$pop1, fst_df$pop2)),
                   value = 0)) %>%
  # Convert to wide matrix
  pivot_wider(names_from = pop2, values_from = value) %>%
  column_to_rownames("pop1") %>%
  as.matrix()

# Sort rows/cols alphabetically
pop_order <- sort(rownames(fst_mat))
fst_mat <- fst_mat[pop_order, pop_order]
print(fst_mat)

# Print matrix
print(fst_mat)


dxy_df <- df %>%
  filter(metric == "Dxy" & !is.na(pop2) & pop2 != ".") %>%
  select(pop1, pop2, value) %>%
  mutate(value = as.numeric(value))

dxy_mat <- dxy_df %>% 
  bind_rows(dxy_df %>% rename(pop1 = pop2, pop2 = pop1)) %>%
  bind_rows(tibble(pop1 = unique(c(dxy_df$pop1, dxy_df$pop2)),
                   pop2 = unique(c(dxy_df$pop1, dxy_df$pop2)),
                   value = 0)) %>%
  pivot_wider(names_from = pop2, values_from = value) %>%
  column_to_rownames("pop1") %>%
  as.matrix()

# Sort rows/cols alphabetically
pop_order <- sort(rownames(dxy_mat))
dxy_mat <- dxy_mat[pop_order, pop_order]
print(dxy_mat)


# ---- 2. Convert to distance object ----
fst_dist <- as.dist(fst_mat)
dxy_dit <- as.dist(dxy_mat)

# ---- 3a. UPGMA clustering ----
fst_hc <- hclust(fst_dist)
dxy_hc <- hclust(dxy_dit)

# Plot quick base R dendrogram
plot(fst_hc, main = "UPGMA clustering (Fst)")
plot(dxy_hc, main = "UPGMA clustering (Dxy)")

# Use the clustering order you computed
fst_ord <- rownames(fst_mat)[fst_hc$order]
dxy_ord <- rownames(dxy_mat)[dxy_hc$order]

pheatmap(
  fst_mat,
  cellwidth=20,
  cellheight=20,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  legend = TRUE, border_color = NA,
  color = rev(paletteer_c("grDevices::Blues", n = 100)),
  display_numbers = F,
  filename = file.path(output_dir, "fst_small_heatmap.png"),
  width = 6,
  height = 6,
)

pheatmap(
  dxy_mat,
  cellwidth=20,
  cellheight=20,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  legend = TRUE, border_color = NA,
  color = rev(paletteer_c("grDevices::Reds", n = 10)),
  display_numbers = F,
  filename = file.path(output_dir, "dxy_small_heatmap.png"),
  width = 6,
  height = 6,
)










