## ==============================================================================
# Batch Processing of Kinematic Data: Signed Angle Calculation and Displacement Normalization
# ==============================================================================
#
# Description:
#   This script performs automated batch processing of motion tracking data
#   stored in CSV format. For each file, it computes:
#
#     1. Signed head angle (degrees) based on a reference vector
#     2. Frame-to-frame displacement of a tracked point (p1)
#     3. Normalized displacement using baseline length (3 × x1)
#
#   The pipeline includes:
#     - Data cleaning (removal of invalid rows and non-numeric values)
#     - Coordinate extraction for predefined landmarks (p1–p12)
#     - Frame-wise kinematic computation
#     - Automatic result export (CSV)
#     - Automatic visualization (PNG plots)
#
# Input:
#   - CSV files containing tracked 2D coordinates (x, y) for multiple points
#   - The first column is assumed to be a timestamp and will be removed
#
# Output:
#   - *_analysis.csv : Processed numerical results
#   - *_angle_normX.png : Plots of head angle and normalized displacement
#
# Key Definitions:
#   - Reference vector: p1(0) - p2(0)
#   - Signed angle: computed via atan2(cross, dot), range (-180°, 180°]
#   - Displacement: Euclidean distance of p1 between consecutive frames
#   - x1: baseline length between p1(0) and p5(0)
#   - Normalized displacement: X(t) / (3 × x1)
#
# Author:
#   Wei Liu, 2026
#
# Notes:
#   - Designed for batch processing of motion capture / tracking datasets
#   - Robust to missing values and irregular file naming
#   - Output filenames are sanitized to avoid system errors
#
# ==============================================================================

## 1) Path Settings
# Set in_dir to your specific folder containing CSV files
in_dir  <- ""   
out_dir <- file.path(in_dir, "batch_out")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## 2) Utility: Clean file names to prevent errors in png() due to special characters
safe_name <- function(x) {
  x <- gsub("%", "pct", x)                 # Replace % with pct
  x <- gsub("[/\\\\?<>\\:*|\"']", "_", x)  # Replace other special characters
  x <- gsub("\\s+", "_", x)                # Replace spaces with underscores
  x
}

## 3) Single File Processing Logic
process_one <- function(fpath) {
  message("Processing: ", basename(fpath))
  
  # Read data
  dat <- tryCatch(read.csv(fpath, check.names = FALSE), error = function(e) NULL)
  if (is.null(dat)) { warning("Failed to read: ", fpath); return(invisible(NULL)) }
  
  # Remove timestamp column (matching your logic)
  if (ncol(dat) < 2) { warning("Insufficient columns: ", fpath); return(invisible(NULL)) }
  data <- dat[, -1, drop = FALSE]
  
  # Extract x,y coordinates for p1~p12 based on current indices
  xy_cols <- c(
    1:2, 4:5, 7:8, 10:11, 13:14, 16:17,
    19:20, 22:23, 25:26, 28:29, 31:32, 34:35
  )
  xy_cols <- xy_cols[xy_cols <= ncol(data)]  # Safety: handle cases with fewer columns
  if (length(xy_cols) < 10) {                
    warning("Insufficient coordinate columns: ", fpath); return(invisible(NULL))
  }
  all_points <- data[, xy_cols, drop = FALSE]
  
  # Force numeric conversion to prevent calculation errors
  for (j in seq_len(ncol(all_points))) {
    all_points[[j]] <- suppressWarnings(as.numeric(all_points[[j]]))
  }
  
  # Remove frames containing NA
  valid_rows <- complete.cases(all_points)
  valid_data <- all_points[valid_rows, , drop = FALSE]
  
  # Remove first two rows sequentially (matching original logic)
  if (nrow(valid_data) >= 2) valid_data <- valid_data[-1, , drop = FALSE]
  if (nrow(valid_data) >= 2) valid_data <- valid_data[-1, , drop = FALSE]
  if (nrow(valid_data) < 2) { warning("Insufficient valid frames: ", fpath); return(invisible(NULL)) }
  
  # Column positions for p1, p2, and p5 (indices: p1=1:2, p2=3:4, p5=9:10)
  needed_idx <- c(1:2, 3:4, 9:10)
  if (max(needed_idx) > ncol(valid_data)) {
    warning("File missing required columns for p1/p2/p5: ", fpath); return(invisible(NULL))
  }
  
  p1_0 <- as.numeric(valid_data[1, 1:2])
  p2_0 <- as.numeric(valid_data[1, 3:4])
  p5_0 <- as.numeric(valid_data[1, 9:10])
  
  # Reference vector: p1 - p2
  reference_vector <- p1_0 - p2_0
  
  # Frame-by-frame calculation (Signed angles + Displacement)
  nF <- nrow(valid_data)
  angles <- numeric(nF - 1)
  displacements <- numeric(nF - 1)
  
  for (i in 2:nF) {
    p1_t <- as.numeric(valid_data[i, 1:2])
    p2_t <- as.numeric(valid_data[i, 3:4])
    
    # Signed angle: atan2(cross, dot) in (-180, 180]
    cur_vec <- p1_t - p2_t
    ref_vec <- reference_vector
    dot   <- sum(ref_vec * cur_vec)
    cross <- ref_vec[1] * cur_vec[2] - ref_vec[2] * cur_vec[1]
    angle_deg <- atan2(cross, dot) * 180 / pi
    
    angles[i - 1] <- angle_deg
    
    # Displacement: p1 relative to previous frame p1
    prev_p1 <- as.numeric(valid_data[i - 1, 1:2])
    displacements[i - 1] <- sqrt(sum((p1_t - prev_p1)^2))
  }
  
  # Normalize displacement: X(t) / (3 * x1), where x1 = ||p1_0 - p5_0||
  x1 <- sqrt(sum((p1_0 - p5_0)^2))
  normalized_displacements <- displacements / (3 * x1)
  
  result <- data.frame(
    Frame = 2:nF,
    Normalized_Displacement_X_t_3_x1 = normalized_displacements,
    Head_Angle_deg = angles
  )
  
  # Generate safe file names
  base <- tools::file_path_sans_ext(basename(fpath))
  base_safe <- safe_name(base)
  
  # Export result table
  csv_out <- file.path(out_dir, paste0(base_safe, "_analysis.csv"))
  write.csv(result, csv_out, row.names = FALSE)
  
  # Plotting and saving (Angle + Normalized Displacement)
  png_out <- file.path(out_dir, paste0(base_safe, "_angle_normX.png"))
  png(png_out, width = 1400, height = 1000, res = 150)
  op <- par(mfrow = c(2,1), mar = c(4,4,2,1))
  plot(result$Frame, result$Head_Angle_deg, type = "b",
       xlab = "Frame", ylab = "Head Angle (degrees)")
  abline(h = 0, lty = 2)
  plot(result$Frame, result$Normalized_Displacement_X_t_3_x1, type = "b",
       xlab = "Frame", ylab = "X(t) / (3 * x1)")
  par(op)
  dev.off()
  
  invisible(result)
}

## 4) Batch Iterate Through CSV Files
files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
out_list <- lapply(files, process_one)

message("Process Complete! Results (CSV + PNG) saved to: ", out_dir)