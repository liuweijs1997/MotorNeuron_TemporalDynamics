# -------------------------------------------------------------------------
# Violin plot of per-cell gene expression (Control vs Training)
#
# This script generates violin plots with overlaid boxplots and jittered
# points to visualize per-cell gene expression levels between experimental
# groups.
#
# The input data should be a per-cell expression matrix (linear scale),
# where rows represent cells and columns include:
#   - gene expression values
#   - a grouping variable (orig.ident)
#
# The script performs:
#   1. Data cleaning
#   2. Reshaping to long format
#   3. Visualization with statistical comparison
#
# Study:
# "Agility Training Enhances Motor Temporal Precision by Reweighting
#  Spinal Phase-Locked Commissural Inhibition"
#
# Dependencies:
#   - dplyr
#   - tidyr
#   - ggplot2
#   - ggpubr
#   - tibble
#
# Author:
#   Wei Liu, 2026
#
# Notes:
#   - Input CSV should be generated from per-cell expression extraction
#   - Statistical comparison uses Wilcoxon rank-sum test
# -------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(tibble)

# Load data (replace with slow or fast dataset)
df <- read.csv(
  "your_file.csv",
  check.names = FALSE
)

# 1) Remove rows with missing group labels
df <- df[!is.na(df$orig.ident), ]

# 2) Remove columns with empty or invalid names
bad_cols <- which(is.na(colnames(df)) | trimws(colnames(df)) == "")
if (length(bad_cols) > 0) {
  df <- df[, -bad_cols, drop = FALSE]
}

# 3) If the first column represents cell/sample IDs, set as row names
if ("Unnamed: 0" %in% colnames(df)) {
  rownames(df) <- df$`Unnamed: 0`
  df$`Unnamed: 0` <- NULL
}

# 4) Specify genes to plot
genes <- c("glra1","glra2","glra3","glra4a","glra4b","glrba","glrbb")

# Ensure all genes exist in the dataset
stopifnot(all(genes %in% colnames(df)))

# 5) Convert to long format for plotting
df_long <- df |>
  rownames_to_column("cell") |>
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "expr"
  ) |>
  mutate(
    gene = factor(gene, levels = genes),
    orig.ident = factor(orig.ident, levels = c("Control", "Training"))
  )

# 6) Define comparison groups
comparisons <- list(c("Control", "Training"))

# 7) Define dodge position for grouped plots
pd <- position_dodge(width = 0.8)

# 8) Generate plot
p <- ggplot(df_long, aes(x = gene, y = expr, fill = orig.ident)) +
  
  # Violin plot: filled without outline
  geom_violin(
    position = pd,
    trim = TRUE,
    scale = "width",
    color = NA
  ) +
  
  # Boxplot: white fill with group-colored outline
  geom_boxplot(
    aes(
      group = interaction(gene, orig.ident),
      color = orig.ident
    ),
    position = pd,
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    size = 0.35
  ) +
  
  # Jittered points: white fill with group-colored border
  geom_point(
    aes(color = orig.ident),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    shape = 21,
    fill = "white",
    stroke = 0.3,
    size = 3,
    alpha = 1
  ) +
  
  # Statistical significance (Wilcoxon test)
  stat_compare_means(
    aes(group = orig.ident),
    method = "wilcox.test",
    label = "p.signif",
    label.y.npc = "top"
  ) +
  
  # Color settings
  scale_fill_manual(values = c(
    "Control" = "#000000",
    "Training" = "#b137af"
  )) +
  
  scale_color_manual(values = c(
    "Control" = "black",
    "Training" = "#b137af"
  )) +
  
  # Axis and layout
  scale_y_continuous(
    limits = c(0, 4),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  
  labs(x = NULL, y = "Expression Level")

# Display plot
p