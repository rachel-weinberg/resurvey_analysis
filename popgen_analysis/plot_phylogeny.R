# plot_phylogeny.R
# Plot SNAPPER MCC tree with clades of posterior < 0.7 collapsed to polytomies,
# tip shape = supercolony, color = timepoint, open = changed / filled = unchanged

pacman::p_load(tidyverse, treeio, ggtree, ape)

# input tree and support threshold
tree_file <- "/Users/rachelweinberg/Documents/UCB/lab_stuff/LhumGenomic/UCEanalysis/phylogeny/snapper_logs/snapper_4run_merged.tree"
min_posterior <- 0.7

# behavior/timepoint pop designations
LH_s1 <- c("LH00", "TWC99")
LH_s2 <- c("LH23", "TWC22")
LH_s2_changed <- c("SB22", "SJ22", "LA22", "LJM22", "LJC22", "PA24")
LSC_unchanged_s1 <- c(
  "KC98",
  "LP04",
  "MB98",
  "MP07",
  "MT98",
  "SAL98",
  "SAU98",
  "SLO98",
  "SM98",
  "UK98",
  "WL03"
)
LSC_changed_s1 <- c("LA98", "LJC00", "LJM98", "PA99", "SB98", "SJ98")
LSC_s2 <- c(
  "EPOW22",
  "UK22",
  "KC23",
  "LP23",
  "MB23",
  "MP23",
  "MT23",
  "SM22",
  "SAL22",
  "SAU22",
  "SLO22",
  "WL22"
)
LS_s1 <- "LS97"
LS_s2 <- "LS23"

beh_tp_df <- bind_rows(
  tibble(pop = LH_s1, beh_tp = "LH_s1"),
  tibble(pop = LH_s2, beh_tp = "LH_s2"),
  tibble(pop = LH_s2_changed, beh_tp = "LH_s2_changed"),
  tibble(pop = LSC_changed_s1, beh_tp = "LSC_changed_s1"),
  tibble(pop = LSC_unchanged_s1, beh_tp = "LSC_unchanged_s1"),
  tibble(pop = LSC_s2, beh_tp = "LSC_s2"),
  tibble(pop = LS_s1, beh_tp = "LS_s1"),
  tibble(pop = LS_s2, beh_tp = "LS_s2")
)

# plot colors
color_beh <- c("#1F78B4", "#FF7F00", "#33A02C", "#6A3D9A", "#B3B3B3")
color_beh_2 <- c("#A6CEE3", "#FDBF6F", "#B2DF8A", "#CAB2D6", "#B15928")
timepoint_colors <- c(
  color_beh[1],
  color_beh_2[1],
  color_beh[2],
  color_beh_2[2],
  color_beh[3],
  color_beh_2[4],
  color_beh[4]
)
# tip style per pop group: shape = supercolony, color = timepoint,
# open = changed / filled = unchanged
tip_style <- tribble(
  ~beh_tp,            ~legend_label,          ~color,              ~shape, ~fill,
  "LSC_unchanged_s1", "LSC s1 unchanged",     timepoint_colors[2], 21,     timepoint_colors[2],
  "LSC_changed_s1",   "LSC s1 changed",       timepoint_colors[2], 21,     "white",
  "LSC_s2",           "LSC s2",               timepoint_colors[1], 21,     timepoint_colors[1],
  "LH_s1",            "LH s1",                timepoint_colors[4], 24,     timepoint_colors[4],
  "LH_s2",            "LH s2 unchanged",      timepoint_colors[3], 24,     timepoint_colors[3],
  "LH_s2_changed",    "LH s2 changed",        timepoint_colors[3], 24,     "white",
  "LS_s1",            "SK s1",                "#a2dd9d",           22,     "#a2dd9d",
  "LS_s2",            "SK s2",                timepoint_colors[5], 22,     timepoint_colors[5],
  "other",            "other",                "grey50",            23,     "grey50"
)
tip_colors <- setNames(tip_style$color, tip_style$legend_label)
tip_shapes <- setNames(tip_style$shape, tip_style$legend_label)
tip_fills <- setNames(tip_style$fill, tip_style$legend_label)

# collapse nodes below support threshold while keeping tip heights;
# each collapsed edge's length is added to its child edges (root-to-tip order)
collapse_low_support <- function(phy, support, threshold) {
  n_tip <- Ntip(phy)
  low_nodes <- n_tip + which(support < threshold)
  low_nodes <- low_nodes[low_nodes != n_tip + 1]
  low_nodes <- intersect(phy$edge[, 2], low_nodes)
  for (nd in low_nodes) {
    parent_edge <- which(phy$edge[, 2] == nd)
    child_edges <- which(phy$edge[, 1] == nd)
    phy$edge.length[child_edges] <- phy$edge.length[child_edges] +
      phy$edge.length[parent_edge]
    phy$edge.length[parent_edge] <- 0
  }
  di2multi(phy, tol = 1e-12)
}

# read tree and attach posterior as node labels;
# treeio mis-parses the snapperCoalescentRate annotations, so strip every annotation except posterior before reading
tree_txt <- readLines(tree_file)
tree_txt <- gsub(
  "\\[&[^]]*?(posterior=[0-9.eE+-]+)[^]]*\\]",
  "[&\\1]",
  tree_txt,
  perl = TRUE
)
tree_txt <- gsub("\\[&(?!posterior=)[^]]*\\]", "", tree_txt, perl = TRUE)
tree_tmp <- tempfile(fileext = ".tree")
writeLines(tree_txt, tree_tmp)

mcc <- read.beast(tree_tmp)
mcc_phylo <- as.phylo(mcc)
node_post <- get.data(mcc) |>
  select(node, posterior) |>
  mutate(node = as.integer(node))
mcc_phylo$node.label <- node_post$posterior[
  match(Ntip(mcc_phylo) + seq_len(Nnode(mcc_phylo)), node_post$node)
]

# collapse clades with posterior < min_posterior
mcc_poly <- collapse_low_support(
  mcc_phylo,
  mcc_phylo$node.label,
  min_posterior
)
mcc_poly$node.label <- round(mcc_poly$node.label, 2)

# tip metadata for shape/color/fill
tip_beh <- tibble(label = mcc_poly$tip.label) |>
  left_join(beh_tp_df, by = c("label" = "pop")) |>
  mutate(beh_tp = replace_na(beh_tp, "other")) |>
  left_join(select(tip_style, beh_tp, legend_label), by = "beh_tp") |>
  mutate(legend_label = factor(legend_label, levels = tip_style$legend_label))

# plot
mcc_plot <- ggtree(mcc_poly) %<+%
  tip_beh +
  geom_tippoint(
    aes(color = legend_label, shape = legend_label, fill = legend_label),
    size = 2.5,
    stroke = 1,
    position = position_nudge(x = 0.002)
  ) +
  geom_tiplab(size = 3, offset = 0.007) +
  geom_nodelab(
    size = 2.5,
    hjust = 1.2,
    vjust = -0.5
  ) +
  scale_color_manual(values = tip_colors, name = "Population group") +
  scale_shape_manual(values = tip_shapes, name = "Population group") +
  scale_fill_manual(values = tip_fills, name = "Population group") +
  theme_tree2() +
  hexpand(0.15) +
  labs(
    title = paste0(
      "SNAPPER MCC tree, clades with posterior < ",
      min_posterior,
      " collapsed"
    )
  )
mcc_plot

ggsave("snapper_mcc_polytomy_beh.svg", mcc_plot, width = 8, height = 10)
