# -------------------------------------------------------------------------
# Per-cell gene expression extraction from Seurat object
#
# This script extracts per-cell expression values for selected genes from
# a Seurat object and converts log-normalized data back to linear scale.
#
# The output is used for downstream visualization (e.g., violin plots)
# and statistical comparison between experimental groups.
#
# Study:
# "Agility Training Enhances Motor Temporal Precision by Reweighting
#  Spinal Phase-Locked Commissural Inhibition"
#
# Dependencies:
#   - Seurat
#
# Author:
#   Wei Liu, 2026
#
# Notes:
#   - Input data should be log1p-normalized (Seurat default: slot = "data")
#   - expm1() is used to convert log1p values back to linear scale
#   - Output is per-cell expression matrix with group labels
#
# -------------------------------------------------------------------------

if (!require(Seurat)) install.packages("Seurat")

# Load Seurat object
MN <- readRDS(file = 'scobj.rds')

# Target genes
genes <- c('glra1','glra2','glra3','glra4a','glra4b','glrba','glrbb')

# ================================
# FAST MOTOR NEURONS
# ================================

obj <- MNsp$`Fast MNs`

# Extract per-cell expression (log1p-normalized values)
mat_log <- GetAssayData(obj, assay = "RNA", slot = "data")[genes, , drop = FALSE]

# Convert log1p values back to linear scale (IMPORTANT)
mat_lin <- expm1(as.matrix(mat_log))

# Convert to data frame (cells as rows)
cell_gene_lin <- as.data.frame(t(mat_lin))

# Add group label
cell_gene_lin$orig.ident <- obj$orig.ident

# Preview
head(cell_gene_lin)

# Save output
write.csv(
  x = cell_gene_lin,
  file = 'fast_gene_expression_per_cell.csv'
)

# ================================
# SLOW MOTOR NEURONS
# ================================

obj <- MNsp$`Slow MNs`

# Extract per-cell expression (log1p-normalized values)
mat_log <- GetAssayData(obj, assay = "RNA", slot = "data")[genes, , drop = FALSE]

# Convert log1p values back to linear scale (IMPORTANT)
mat_lin <- expm1(as.matrix(mat_log))

# Convert to data frame (cells as rows)
cell_gene_lin <- as.data.frame(t(mat_lin))

# Add group label
cell_gene_lin$orig.ident <- obj$orig.ident

# Preview
head(cell_gene_lin)

# Save output
write.csv(
  x = cell_gene_lin,
  file = 'slow_gene_expression_per_cell.csv'
)