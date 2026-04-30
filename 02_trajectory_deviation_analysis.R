# ==============================================================================
# Batch Polar Analysis of Kinematic Data with Height Projection and Summary Statistics
# ==============================================================================
#
# Description:
#   This script performs batch processing of kinematic analysis outputs
#   (e.g., *_analysis.csv files) to generate polar coordinate visualizations
#   and compute derived height-based metrics.
#
#   For each input file, the pipeline:
#     - Reads time, displacement (distance), and angle data
#     - Converts angular values from degrees to radians
#     - Projects displacement onto the vertical axis (height component)
#     - Computes both signed and absolute height values
#     - Normalizes height into percentage scale (×100)
#     - Generates polar "spoke" plots for visualization
#     - Aggregates summary statistics across all files
#
# Input:
#   - CSV files matching pattern "*_analysis.csv"
#   - Each file must contain at least three columns:
#       time, distance, angle (degrees)
#
# Output:
#   - SVG plots: Polar spoke visualization for each file
#   - Summary CSV:
#       mean_height_abs,
#       mean_height_signed,
#       mean_height_abs_pct,
#       distance and angle ranges
#
# Key Definitions:
#   - angle_rad = angle_deg × π / 180
#   - height_signed = distance × sin(angle_rad)
#   - height_abs = |height_signed|
#   - height_abs_pct = height_abs × 100
#
# Author:
#   Wei Liu, 2026
#
# Notes:
#   - Polar plots use angle as angular coordinate and displacement as radius
#   - Output filenames are automatically generated per input file
#   - Designed for batch processing of motion/behavior tracking data
#
# ==============================================================================


library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(tools)
library(svglite)

dir_in  <- ""
pattern <- "_analysis\\.csv$"
save_plots <- TRUE
out_plot_dir <- file.path(dir_in, "polar_plots_svg")
out_summary  <- file.path(dir_in, "mean_height_abs_pct_summary.csv")

if (save_plots && !dir.exists(out_plot_dir)) dir.create(out_plot_dir, recursive = TRUE)

deg2rad <- function(x) x * pi / 180

plot_spokes <- function(data, value_col = "distance",
                        title = "Polar spokes (r in [0, 0.1])") {
  
  data$distance_plot <- pmin(pmax(data[[value_col]], 0), 0.1)
  
  ggplot(data, aes(x = angle_deg)) +
    geom_segment(aes(xend = angle_deg, y = 0, yend = distance_plot), linewidth = 0.4) +
    geom_point(aes(y = distance_plot), size = 1.2) +
    coord_polar(theta = "x", start = -pi/2, direction = 1, clip = "off") +
    scale_x_continuous(breaks = c(-180, 0, 180), limits = c(-180, 180), minor_breaks = NULL) +
    scale_y_continuous(limits = c(0, 0.1)) +
    labs(x = "Angle (deg, -90…90)", y = NULL, title = title) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

process_one <- function(csv_file, save_plot = TRUE) {
  df <- suppressMessages(read_csv(csv_file, col_types = cols(.default = col_double())))
  if (ncol(df) < 3) stop("Less than 3 columns: ", csv_file)
  df <- df[, 1:3]
  names(df) <- c("time", "distance", "angle_deg")
  
  df <- df %>% filter(!is.na(distance), !is.na(angle_deg))
  if (nrow(df) == 0) stop("No valid rows: ", csv_file)
  
  df <- df %>%
    mutate(angle_rad = deg2rad(angle_deg),
           height_signed = distance * sin(angle_rad),
           height_abs = abs(height_signed),
           height_abs_pct = height_abs * 100)
  
  if (save_plot) {
    p_title <- paste0("Polar spokes per frame (raw distance, r∈[0,0.1])\n", basename(csv_file))
    p <- plot_spokes(df, "distance", p_title)
    svg_out <- file.path(out_plot_dir,
                         paste0(file_path_sans_ext(basename(csv_file)), "_polar.svg"))
    ggsave(svg_out, p, width = 7, height = 7, dpi = 200, device = svglite)
  }
  
  tibble(
    file = basename(csv_file),
    path = csv_file,
    n_rows = nrow(df),
    mean_height_abs = mean(df$height_abs, na.rm = TRUE),
    mean_height_signed = mean(df$height_signed, na.rm = TRUE),
    mean_height_abs_pct = mean(df$height_abs_pct, na.rm = TRUE),
    distance_min = min(df$distance, na.rm = TRUE),
    distance_max = max(df$distance, na.rm = TRUE),
    angle_min = min(df$angle_deg, na.rm = TRUE),
    angle_max = max(df$angle_deg, na.rm = TRUE)
  )
}

files <- list.files(dir_in, pattern = pattern, full.names = TRUE)
if (length(files) == 0) stop("No matching files found in directory: ", dir_in)

results <- purrr::map_dfr(
  files,
  ~ tryCatch(process_one(.x, save_plot = save_plots),
             error = function(e) {
               warning(sprintf("Skipped (error): %s -> %s", .x, e$message))
               tibble(
                 file = basename(.x), path = .x, n_rows = NA_integer_,
                 mean_height_abs = NA_real_, mean_height_signed = NA_real_,
                 mean_height_abs_pct = NA_real_, distance_min = NA_real_,
                 distance_max = NA_real_, angle_min = NA_real_, angle_max = NA_real_
               )
             })
)

readr::write_csv(results, out_summary)
cat("Done! SVG plots saved in: ", out_plot_dir, "\n")
cat("Summary table: ", out_summary, "\n")