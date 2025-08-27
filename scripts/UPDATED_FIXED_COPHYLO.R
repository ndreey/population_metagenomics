library(tidyverse)
library(ape)
library(phytools)
library(svglite)

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

# Load in trees
tree_fly <- read.tree("data/001_phylo/complete_single_diptera_buscos.astral.tre")
tree_stam <- read.tree("data/001_phylo/stam.filt.bial.passed.no-sc.no-outlier.treefile")


# Remove branch lengths
tree_fly$edge.length  <- NULL
tree_stam$edge.length <- NULL

#plot(tree_fly)
#plot(tree_stam)


# Add root to same sample as fly tree
#tree_stam <- root(tree_stam, outgroup = "CHSK_P12002_101", resolve.root = TRUE)

# Extract tip labels
tips_fly  <- tree_fly$tip.label
tips_stam <- tree_stam$tip.label

# Extract sample ID from stam (everything after first underscore)
sample_ids_stam <- str_extract(tips_stam, "P\\d+_\\d+")

# Extract population prefix (everything before that)
pop_prefix_stam <- str_extract(tips_stam, "^[A-Z]+(?=_P\\d+_\\d+)")

# Build a lookup table: SAMPLEID -> POP
lookup <- tibble(
  sample_id = sample_ids_stam,
  pop       = pop_prefix_stam
) %>% distinct()

# Add population info to tree_fly labels (if available)
tips_fly_new <- tips_fly %>%
  map_chr(~ {
    if (.x %in% lookup$sample_id) {
      paste0(lookup$pop[lookup$sample_id == .x], "_", .x)
    } else {
      # keep as is if no match
      .x
    }
  })

# Update the tree
tree_fly$tip.label <- tips_fly_new

# tips in common
common_tips <- intersect(tree_fly$tip.label, tree_stam$tip.label)

# prune tree_fly to only those common tips
tree_fly_pruned <- keep.tip(tree_fly, common_tips)


# Add "symb." prefix to stam tips
stam_tips <- paste0("symb.", tree_stam$tip.label)

# Rename stam tree tips
tree_stam$tip.label <- stam_tips

# Add symb to association
stam_tips_asso <- paste0("symb.", common_tips)

# Association
shared <- as.data.frame(cbind("Stammerula" = stam_tips_asso,
                              "T.conura" =common_tips))

# Visualize both trees side by side
co <- cophylo(tree_stam, tree_fly_pruned, assoc = shared, rotate = T)

# Lets add some meta
meta <- shared %>%
  mutate(POP = str_extract(T.conura, "^[^_]+")) %>%
  left_join(metadata, by = "POP")


# Order of the nodes
left_tips  <- co$trees[[1]]$tip.label      # order the left tree is drawn
right_tips <- co$trees[[2]]$tip.label      # order the right tree is drawn

# map POP color & shape to each tip (left tree labels match meta$T.conura)
right_cols <- setNames(meta$color, meta$T.conura)[right_tips]
right_pch  <- setNames(meta$shape, meta$T.conura)[right_tips]

# right tree labels are "symb." prefixed in meta$Stammerula
left_cols <- setNames(meta$color, meta$Stammerula)[left_tips]
left_pch  <- setNames(meta$shape, meta$Stammerula)[left_tips]



## PNG (high-res)
png("00_thesis-fig/lalala_cophylo_plot.png", width=2250, height=2625, 
    units = "px", res = 300)
par(lend = 3, mar = c(1, 1, 2, 1))
plot(co,
     fsize = 0.45,
     link.col = make.transparent(meta$color, 0.8),
     link.lwd = 3,
     link.type = "curved",
     link.lty = "solid")

## Symbols (shape + color)
tiplabels.cophylo(which="left",  pch=left_pch,  bg=left_cols,  col="black", cex=1)
tiplabels.cophylo(which="right", pch=right_pch, bg=right_cols, col="black", cex=1)
## Branchsupport
#nodelabels.cophylo(which="left", frame="none", cex=0.4)
#nodelabels.cophylo(which="right", frame="none", cex=0.4)
dev.off()

# SAve as SVG
svglite::svglite("00_thesis-fig/lalala_cophylo_plot.svg", width = 7.5, height = 8.75)  # inches
par(lend = 3, mar = c(1, 1, 2, 1))

plot(co,
     fsize = 0.45,
     link.col = make.transparent(meta$color, 0.8),
     link.lwd = 3, link.type="curved", link.lty="solid")

tiplabels.cophylo(which="left",  pch=left_pch,  bg=left_cols,  col="black", cex=1)
tiplabels.cophylo(which="right", pch=right_pch, bg=right_cols, col="black", cex=1)
#nodelabels.cophylo(which="left",  frame="none", cex=0.4)
#nodelabels.cophylo(which="right", frame="none", cex=0.4)

dev.off()

