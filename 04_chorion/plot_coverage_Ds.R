# --- Coverage plots for D. simulans ---
# Reads 1kb binned coverage files, normalizes by per-chromosome mean,
# and produces one plot per sample plus a summary plot per chromosome.
# NOTE: Normalization is per-chromosome. If a duplication spans multiple
# chromosomes or you want genome-wide normalization, compute the global
# mean across all chromosomes first and pass it in as a scaling factor.

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No chromosome specified. Usage: Rscript plot_coverage_Ds.R <chromosome_name>", call. = FALSE)
}
chromosome_name <- args[1]
cat("--- Plotting coverage for chromosome:", chromosome_name, "---\n")

coverage_dir <- "process/Ds_coverage/"
plot_dir     <- "process/Ds_coverage/"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("Dsim_1-L1F", "Dsim_1-L2F", "Dsim_1-L3F",
             "Dsim_2-L1F", "Dsim_2-L2F", "Dsim_2-L3F")

# Color by selection class
sample_colors <- c(
  "Dsim_1-L1F" = "grey50",
  "Dsim_1-L2F" = "steelblue",
  "Dsim_1-L3F" = "sienna",
  "Dsim_2-L1F" = "grey20",
  "Dsim_2-L2F" = "dodgerblue",
  "Dsim_2-L3F" = "darkorange"
)

# --- 1. READ AND NORMALIZE ALL SAMPLES ---
all_data <- list()

for (sample in samples) {
  infile <- paste0(coverage_dir, sample, "_", chromosome_name, "_1kb.txt")

  if (!file.exists(infile)) {
    cat("WARNING: File not found, skipping:", infile, "\n")
    next
  }

  df <- read.table(infile, header = TRUE, sep = "\t")

  # Normalize by per-chromosome mean (duplication shows up as ~2)
  chr_mean <- mean(df$MEAN_DEPTH, na.rm = TRUE)
  if (chr_mean == 0) {
    cat("WARNING: Mean depth is 0 for", sample, "on", chromosome_name, "- skipping\n")
    next
  }
  df$NORM_DEPTH <- df$MEAN_DEPTH / chr_mean
  df$SAMPLE     <- sample

  all_data[[sample]] <- df

  # --- Per-sample plot (capped at 5) ---
  p <- ggplot(df, aes(x = WIN_START / 1e6, y = NORM_DEPTH)) +
    geom_line(color = sample_colors[sample], linewidth = 0.4) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = 2, linetype = "dotted", color = "red",   linewidth = 0.4) +
    scale_x_continuous(name = paste0(chromosome_name, " position (Mb)")) +
    scale_y_continuous(name = "Normalized coverage", limits = c(0, 5)) +
    ggtitle(paste(sample, "-", chromosome_name)) +
    theme_classic(base_size = 11)

  plot_file <- paste0(plot_dir, sample, "_", chromosome_name, "_coverage.png")
  ggsave(plot_file, p, width = 10, height = 3, dpi = 150)
  cat("Saved per-sample plot:", plot_file, "\n")

  # --- Per-sample zoomed plot of high-coverage regions ---
  spike_threshold <- 5
  spike_windows <- df[df$NORM_DEPTH > spike_threshold, ]
  if (nrow(spike_windows) > 0) {
    for (spike_start in unique(floor(spike_windows$WIN_START / 1e6))) {
      zoom_min <- spike_start - 0.5
      zoom_max <- spike_start + 0.5
      zoom_df  <- df[df$WIN_START / 1e6 >= zoom_min & df$WIN_START / 1e6 <= zoom_max, ]
      p_zoom <- ggplot(zoom_df, aes(x = WIN_START / 1e6, y = NORM_DEPTH)) +
        geom_line(color = sample_colors[sample], linewidth = 0.4) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
        geom_hline(yintercept = 2, linetype = "dotted", color = "red",   linewidth = 0.4) +
        scale_x_continuous(name = paste0(chromosome_name, " position (Mb)")) +
        scale_y_continuous(name = "Normalized coverage", limits = c(0, NA)) +
        ggtitle(paste(sample, "-", chromosome_name, "- zoom", zoom_min, "-", zoom_max, "Mb")) +
        theme_classic(base_size = 11)
      zoom_file <- paste0(plot_dir, sample, "_", chromosome_name, "_zoom_", spike_start, "Mb_coverage.png")
      ggsave(zoom_file, p_zoom, width = 10, height = 3, dpi = 150)
      cat("Saved zoomed plot:", zoom_file, "\n")
    }
  }
}

# --- 2. SUMMARY PLOT (all samples overlaid) ---
if (length(all_data) == 0) {
  stop("No data loaded. Check coverage_dir and file names.")
}

combined_df <- do.call(rbind, all_data)

p_summary <- ggplot(combined_df, aes(x = WIN_START / 1e6, y = NORM_DEPTH,
                                      color = SAMPLE, group = SAMPLE)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 2, linetype = "dotted", color = "red",   linewidth = 0.4) +
  scale_color_manual(values = sample_colors, name = "Sample") +
  scale_x_continuous(name = paste0(chromosome_name, " position (Mb)")) +
  scale_y_continuous(name = "Normalized coverage", limits = c(0, 5)) +
  ggtitle(paste("All samples -", chromosome_name, "(capped at 5)")) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

summary_file <- paste0(plot_dir, "summary_", chromosome_name, "_coverage.png")
ggsave(summary_file, p_summary, width = 12, height = 4, dpi = 150)
cat("Saved summary plot:", summary_file, "\n")

# --- Summary zoomed plot of high-coverage regions ---
spike_threshold <- 5
spike_windows <- combined_df[combined_df$NORM_DEPTH > spike_threshold, ]
if (nrow(spike_windows) > 0) {
  for (spike_start in unique(floor(spike_windows$WIN_START / 1e6))) {
    zoom_min <- spike_start - 0.5
    zoom_max <- spike_start + 0.5
    zoom_df  <- combined_df[combined_df$WIN_START / 1e6 >= zoom_min & combined_df$WIN_START / 1e6 <= zoom_max, ]
    p_zoom <- ggplot(zoom_df, aes(x = WIN_START / 1e6, y = NORM_DEPTH,
                                   color = SAMPLE, group = SAMPLE)) +
      geom_line(linewidth = 0.4, alpha = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
      geom_hline(yintercept = 2, linetype = "dotted", color = "red",   linewidth = 0.4) +
      scale_color_manual(values = sample_colors, name = "Sample") +
      scale_x_continuous(name = paste0(chromosome_name, " position (Mb)")) +
      scale_y_continuous(name = "Normalized coverage", limits = c(0, NA)) +
      ggtitle(paste("All samples -", chromosome_name, "- zoom", zoom_min, "-", zoom_max, "Mb")) +
      theme_classic(base_size = 11) +
      theme(legend.position = "right")
    zoom_file <- paste0(plot_dir, "summary_", chromosome_name, "_zoom_", spike_start, "Mb_coverage.png")
    ggsave(zoom_file, p_zoom, width = 12, height = 4, dpi = 150)
    cat("Saved zoomed summary plot:", zoom_file, "\n")
  }
}

# --- 3. PEAK REGION PLOTS FOR chr3L ---
if (chromosome_name == "chr3L") {
  peak_regions <- list(
    list(name = "fullRegion_7.85to8.83Mb",  start = 7851416, end = 8834090),
    list(name = "peakRegion_8.1to8.2Mb",    start = 8100000, end = 8200000)
  )
  for (region in peak_regions) {
    zoom_df <- combined_df[combined_df$WIN_START >= region$start &
                           combined_df$WIN_END   <= region$end, ]
    p_peak <- ggplot(zoom_df, aes(x = WIN_START / 1e6, y = NORM_DEPTH,
                                   color = SAMPLE, group = SAMPLE)) +
      geom_line(linewidth = 0.4, alpha = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
      geom_hline(yintercept = 2, linetype = "dotted", color = "red",   linewidth = 0.4) +
      scale_color_manual(values = sample_colors, name = "Sample") +
      scale_x_continuous(name = "chr3L position (Mb)") +
      scale_y_continuous(name = "Normalized coverage", limits = c(0, 5)) +
      ggtitle(paste("All samples - chr3L -", region$name)) +
      theme_classic(base_size = 11) +
      theme(legend.position = "right")
    peak_file <- paste0(plot_dir, "summary_chr3L_", region$name, "_coverage.png")
    ggsave(peak_file, p_peak, width = 12, height = 4, dpi = 150)
    cat("Saved peak region plot:", peak_file, "\n")
  }
}

cat("--- Coverage plotting complete for chromosome:", chromosome_name, "---\n")

# --- 4. ZIP ALL PLOTS IF ALL CHROMOSOMES ARE DONE ---
all_chrs <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
summary_files <- paste0(plot_dir, "summary_", all_chrs, "_coverage.png")

if (all(file.exists(summary_files))) {
  zip_file <- paste0(plot_dir, "Ds_coverage_plots.zip")
  all_plots <- list.files(plot_dir, pattern = "\\.png$", full.names = TRUE)
  zip(zip_file, files = all_plots, flags = "-j")
  cat("All chromosomes complete. Plots zipped to:", zip_file, "\n")
}
