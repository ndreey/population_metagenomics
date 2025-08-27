library(tidyverse)
library(ape)
library(ggtree)

# --- metadata & palettes (as in your script) ---
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
  HP = c("CH" = "#4A1486","CO" = "#00441B"),
  REGION = c("west" = "wheat","east" = "wheat4","scot" = "steelblue1"),
  regHP = c("westCH"="orchid1","eastCH"="orchid4","westCO"="palegreen1","eastCO"="palegreen4","scotHP"="skyblue3"),
  POP = c(
    "CHES"="#BCBDDC","CHFI"="#6A51A3","CHSC"="cadetblue","CHSK"="#9E9AC8",
    "CHST"="#4A1486","COES"="#C7E9C0","COGE"="#00441B","COLI"="#238B45",
    "COSK"="#74C476","CPSC"="cadetblue1")
)

p_shape <- list(
  POP = c(
    "CHES"=23,"CHFI"=23,"CHSC"=20,"CHSK"=21,"CHST"=21,
    "COES"=23,"COGE"=21,"COLI"=23,"COSK"=21,"CPSC"=20)
)

# --- load tree ---
tree <- read.tree("data/001_phylo/stam.filt.bial.passed.no-sc.no-outlier.treefile")

# --- extract POP from tip labels and join metadata ---
tip_annot <- tibble(label = tree$tip.label) |>
  mutate(POP = sub("_.*$", "", label)) |>
  left_join(metadata, by = "POP")

p <- ggtree(tree, layout = "fan", open.angle = 60, branch.length = "none") %<+% tip_annot +
  # tip symbols
  geom_tippoint(aes(shape = POP, fill = POP), size = 3, stroke = 0.7) +
  scale_fill_manual(values = p_colors$POP, name = "Population") +
  scale_shape_manual(values = p_shape$POP, name = "Population") +
  
  # POP text at tips
  geom_tiplab2(
    aes(label = POP),
    color = "black",
    align = FALSE,
    offset = 3,
    hjust = 0.5,
    size = 3,
    show.legend = FALSE
  ) +
  
  # bootstrap as colored dots (discrete bins)
  geom_point2(
    aes(
      subset = !isTip & !is.na(as.numeric(label)),
      color = cut(
        as.numeric(label),
        breaks = c(-Inf, 49, 69, 89, 99, 100),
        labels = c("<50", "50–69", "70–89", "90–99", "100")
      )
    ),
    size = 1
  ) +
  scale_color_manual(
    name = "Bootstrap",
    values = c(
      "<50"   = "#1B1919",
      "50–69" = "#E18727",
      "70–89" = "#0072B5",
      "90–99" = "#20854E",
      "100"   = "#BC3C29"
    )
  ) +
  
  theme_tree() +
  theme(legend.position = "bottom")

p



# save plot as png
ggsave("00_thesis-fig/STAM_PHYLO_CIRC_label.png", p, width = 7, height = 7, dpi = 300)
# Save plot as svg
ggsave("00_thesis-fig/STAM_PHYLO_CIRC._label.svg", p, width = , height = 7)



p2 <- ggtree(tree, layout = "rectangular", size = 1) %<+% tip_annot +
  # tip symbols
  geom_tippoint(aes(shape = POP, fill = POP), size = 3, stroke = 0.7) +
  scale_fill_manual(values = p_colors$POP, name = "Population") +
  scale_shape_manual(values = p_shape$POP, name = "Population") +
  
  # POP text at tips
  geom_tiplab(
    aes(label = POP),
    color = "black",
    align = F,
    hjust = -0.5,
    size = 3,
    show.legend = FALSE
  ) +
  
  # bootstrap as colored dots (discrete bins)
  geom_point2(
    aes(
      subset = !isTip & !is.na(as.numeric(label)),
      color = cut(
        as.numeric(label),
        breaks = c(-Inf, 49, 69, 89, 99, 100),
        labels = c("<50", "50–69", "70–89", "90–99", "100")
      )
    ),
    size = 1
  ) +
  scale_color_manual(
    name = "Bootstrap",
    values = c(
      "<50"   = "#1B1919",
      "50–69" = "#E18727",
      "70–89" = "#0072B5",
      "90–99" = "#20854E",
      "100"   = "#BC3C29"
    )
  ) +
  
  theme_tree2() +
  theme(legend.position = "right")

p2

# save plot as png
ggsave("00_thesis-fig/STAM_PHYLO_RECT_label.png", p2, width = 7, height = 7, dpi = 300)
# Save plot as svg
ggsave("00_thesis-fig/STAM_PHYLO_RECT_label.svg", p2, width = 7, height = 7)



p3 <- ggtree(tree, layout = "rectangular", size = 0.6) %<+% tip_annot +
  # tip symbols
  geom_tippoint(aes(shape = POP, fill = POP), size = 3, stroke = 0.7) +
  scale_fill_manual(values = p_colors$POP, name = "Population") +
  scale_shape_manual(values = p_shape$POP, name = "Population") +
  
  # POP text at tips
  geom_tiplab(
    aes(label = POP),
    color = "black",
    align = FALSE,
    hjust = -0.3,
    size = 3,
    show.legend = FALSE
  ) +
  
  # bootstrap values as text
  geom_text2(
    aes(
      subset = !isTip & !is.na(as.numeric(label)),
      label = label        # directly show bootstrap values
    ),
    size = 2.5,             # adjust text size
    hjust = 1.5,           # tweak positioning if overlapping
    vjust = -0.5,
    color = "black"
  ) +
  
  theme_tree2() +
  theme(legend.position = "right")
p3

# save plot as png
ggsave("00_thesis-fig/STAM_PHYLO_RECT_label_bs.png", p3, width = 10, height = 13, dpi = 300)
# Save plot as svg
ggsave("00_thesis-fig/STAM_PHYLO_RECT_label_bs.svg", p3, width = 10, height = 13)

