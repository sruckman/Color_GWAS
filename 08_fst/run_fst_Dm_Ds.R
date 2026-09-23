#!/usr/bin/env Rscript
# run_fst_Dm_Ds.R
#
# Fst between the two CONTROL pools in each species, windowed along the genome.
#
# Comparisons:
#   D. mel: Toronto (HOULE_L1F) vs Miami (HOUSTON_L1F)
#   D. sim: Tallahassee-1 (Dsim_1-L1F) vs Tallahassee-2 (Dsim_2-L1F)
#
# Prediction:
#   Toronto vs Miami: substantial Fst with hotspots, especially over the
#     cosmopolitan inversions on 2L, 3L, and 3R that Miami carries but Toronto does not.
#   Tallahassee-1 vs Tallahassee-2: low Fst, essentially flat, because they are
#     temporal samples of the same locality drawing on the same standing variation.
#
# Estimator: Hudson's Fst (Bhatia et al. 2013), ratio-of-sums form for windowing.
#   p1 = ALT1 / (REF1 + ALT1);   n1 = REF1 + ALT1
#   p2 = ALT2 / (REF2 + ALT2);   n2 = REF2 + ALT2
#   num   = (p1-p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
#   denom = p1*(1-p2) + p2*(1-p1)
#   window Fst = sum(num) / sum(denom)   [per-SNP Fst = num/denom]
#
# RefAlt source: tries original non-downsampled files first; falls back to
#   downsampled Dark_vs_Control files if the originals are not present.
#   For D. simulans, the original RefAlt files use NCBI chromosome names
#   (NC_052520.2 etc.), so the downsampled chr-named files are preferred.
#
# Output: process/down_sample/bam_ds/plots_censored/fst/

suppressMessages(library(tidyverse))

# ============================================================
# PATHS
# ============================================================
dm_orig  <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/Dm_color"
dm_ds    <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/Dark_vs_Control/Dm_color_sep"
ds_ds    <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds/Dark_vs_Control/Ds_color_sep"

out_dir  <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/fst"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chroms   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
MIN_DEPTH <- 20L    # minimum reads per pool at a SNP
WIN_10    <- 10000L
WIN_50    <- 50000L
MIN_SNPS  <- 20L    # minimum SNPs per window

# Chorion regions to exclude (same windows as the Manhattan script)
chorion_mel <- tribble(
  ~chr,    ~start,    ~end,
  "chr3L", 8500000L,  9000000L,
  "chrX",  8400000L,  8600000L
)
chorion_sim <- tribble(
  ~chr,    ~start,    ~end,
  "chr3L", 8342753L,  8834090L,
  "chrX",  7800000L,  8400000L   # Confirmed by coverage scan. Peak at 8,208,001 (up to 68x).
)

# D. mel inversion windows (all 4 cosmopolitan inversions segregating in Miami)
# Coordinates: dm6, matching plot_bam_ds_manhattan.R
inversions_mel <- tribble(
  ~name,       ~chr,    ~start,     ~end,
  "In(2L)t",   "chr2L",  2225744L,  13154180L,
  "In(3L)P",   "chr3L",  3173046L,  16301941L,
  "In(3R)K",   "chr3R",  7576289L,  22422742L,
  "In(3R)P",   "chr3R", 12257931L,  21082440L
)

# ============================================================
# FUNCTIONS
# ============================================================

read_refalt <- function(dir, chr) {
  f <- file.path(dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(f)) return(NULL)
  read.table(f, header = TRUE, sep = "", check.names = FALSE) |>
    as_tibble()
}

hudson_fst <- function(ref1, alt1, ref2, alt2) {
  n1 <- ref1 + alt1
  n2 <- ref2 + alt2
  p1 <- alt1 / n1
  p2 <- alt2 / n2
  num   <- (p1 - p2)^2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
  denom <- p1 * (1 - p2) + p2 * (1 - p1)
  list(num = num, denom = denom)
}

# Read one chromosome from a RefAlt file and compute per-SNP Hudson numerator/denominator.
# col1/col2 are the sample codes (e.g. "HOULE_L1F", "HOUSTON_L1F").
load_chr_fst <- function(refalt_dir, chr, col1, col2,
                          chorion_tbl, min_depth = MIN_DEPTH) {
  d <- read_refalt(refalt_dir, chr)
  if (is.null(d)) { cat(sprintf("  [missing] RefAlt.%s.txt\n", chr)); return(NULL) }

  ref1 <- d[[paste0("REF_", col1)]]
  alt1 <- d[[paste0("ALT_", col1)]]
  ref2 <- d[[paste0("REF_", col2)]]
  alt2 <- d[[paste0("ALT_", col2)]]

  if (is.null(ref1) || is.null(ref2))
    stop(sprintf("Columns not found for %s or %s in %s", col1, col2, chr))

  fst_vals <- hudson_fst(ref1, alt1, ref2, alt2)

  out <- tibble(
    chr   = chr,
    pos   = d$POS,
    n1    = ref1 + alt1,
    n2    = ref2 + alt2,
    fst_num   = fst_vals$num,
    fst_denom = fst_vals$denom
  ) |>
    filter(n1 >= min_depth, n2 >= min_depth, is.finite(fst_num), fst_denom > 0)

  # Remove chorion regions
  for (i in seq_len(nrow(chorion_tbl))) {
    out <- out |> filter(!(chr == chorion_tbl$chr[i] &
                             pos >= chorion_tbl$start[i] &
                             pos <= chorion_tbl$end[i]))
  }
  out
}

window_fst <- function(dat, win_size, min_snps = MIN_SNPS) {
  dat |>
    mutate(win = floor((pos - 1) / win_size) * win_size + 1L) |>
    group_by(chr, win) |>
    summarise(
      n_snps   = n(),
      fst      = sum(fst_num) / sum(fst_denom),
      .groups  = "drop"
    ) |>
    filter(n_snps >= min_snps, is.finite(fst)) |>
    mutate(win_mid = win + win_size / 2)
}

hotspots <- function(wins, top_pct = 0.01) {
  thresh <- quantile(wins$fst, 1 - top_pct, na.rm = TRUE)
  wins |> filter(fst >= thresh) |> arrange(desc(fst))
}

# ============================================================
# LOAD DATA
# ============================================================

# --- D. mel: try original RefAlt, fall back to downsampled ---
test_f <- file.path(dm_orig, "RefAlt.chrX.txt")
dm_dir <- if (file.exists(test_f)) {
  cat("Using original (non-downsampled) D. mel RefAlt files.\n")
  dm_orig
} else {
  cat("Original D. mel RefAlt not found; using downsampled files.\n")
  dm_ds
}

cat("\nLoading D. mel (Toronto vs Miami) ...\n")
dm_snps <- map_dfr(chroms, ~ load_chr_fst(dm_dir, .x,
                     "HOULE_L1F", "HOUSTON_L1F", chorion_mel))
cat(sprintf("  %d SNPs loaded after depth and chorion filters.\n", nrow(dm_snps)))

# --- D. sim: use downsampled (original uses NCBI chromosome names) ---
cat("\nLoading D. sim (Tallahassee-1 vs Tallahassee-2) ...\n")
ds_snps <- map_dfr(chroms, ~ load_chr_fst(ds_ds, .x,
                     "Dsim_1-L1F", "Dsim_2-L1F", chorion_sim))
cat(sprintf("  %d SNPs loaded after depth and chorion filters.\n", nrow(ds_snps)))

# ============================================================
# WINDOW FST
# ============================================================
dm_10 <- window_fst(dm_snps, WIN_10)
dm_50 <- window_fst(dm_snps, WIN_50)
ds_10 <- window_fst(ds_snps, WIN_10)
ds_50 <- window_fst(ds_snps, WIN_50)

write_csv(dm_50, file.path(out_dir, "fst_Toronto_vs_Miami_50kb.csv"))
write_csv(ds_50, file.path(out_dir, "fst_Tallahassee1_vs_Tallahassee2_50kb.csv"))
write_csv(dm_10, file.path(out_dir, "fst_Toronto_vs_Miami_10kb.csv"))
write_csv(ds_10, file.path(out_dir, "fst_Tallahassee1_vs_Tallahassee2_10kb.csv"))

# ============================================================
# SUMMARY STATISTICS
# ============================================================
cat("\n########## Genome-wide windowed Fst (50 kb windows) ##########\n")
bind_rows(
  dm_50 |> summarise(comparison = "Toronto vs Miami (D. mel)",
                     n_windows = n(), mean_fst = mean(fst), median_fst = median(fst),
                     pct99 = quantile(fst, 0.99), pct999 = quantile(fst, 0.999)),
  ds_50 |> summarise(comparison = "Tallahassee-1 vs Tallahassee-2 (D. sim)",
                     n_windows = n(), mean_fst = mean(fst), median_fst = median(fst),
                     pct99 = quantile(fst, 0.99), pct999 = quantile(fst, 0.999))
) |> as.data.frame() |> print()

cat("\n########## Top 1% Fst hotspots — Toronto vs Miami (D. mel, 50 kb) ##########\n")
dm_hot1 <- hotspots(dm_50, 0.01)
dm_hot1 |> as.data.frame() |> head(20) |> print()
write_csv(dm_hot1, file.path(out_dir, "hotspots_Toronto_vs_Miami_top1pct.csv"))

cat("\n########## Top 1% Fst hotspots — Tallahassee (D. sim, 50 kb) ##########\n")
ds_hot1 <- hotspots(ds_50, 0.01)
ds_hot1 |> as.data.frame() |> head(20) |> print()
write_csv(ds_hot1, file.path(out_dir, "hotspots_Tallahassee_top1pct.csv"))

# Flag D. mel hotspots that overlap the cosmopolitan inversions
cat("\n########## D. mel hotspots overlapping cosmopolitan inversions ##########\n")
inv_hits <- dm_hot1 |>
  rowwise() |>
  filter(any(
    chr == inversions_mel$chr &
    win_mid >= inversions_mel$start &
    win_mid <= inversions_mel$end
  )) |>
  ungroup()
if (nrow(inv_hits) > 0) {
  cat(sprintf("  %d of %d top-1%% hotspot windows fall inside an inversion.\n",
              nrow(inv_hits), nrow(dm_hot1)))
  inv_hits |> as.data.frame() |> print()
} else {
  cat("  No top-1% hotspot windows fall inside the listed inversions.\n")
}

# ============================================================
# PLOTS — single horizontal cumulative-position panel
# ============================================================
chrom_order <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
GAP_MB      <- 2e6L   # visual gap between chromosomes

plot_fst <- function(wins, title, win_size = WIN_50, inv_tbl = NULL) {
  # Cumulative chromosome offsets from data
  chr_sizes <- wins |>
    group_by(chr) |>
    summarise(size = max(win_mid) + win_size, .groups = "drop") |>
    mutate(chr = factor(chr, levels = chrom_order)) |>
    arrange(chr) |>
    mutate(offset = cumsum(lag(size, default = 0L)) +
                    (as.integer(chr) - 1L) * GAP_MB)

  off <- setNames(chr_sizes$offset, as.character(chr_sizes$chr))

  d <- wins |>
    filter(chr %in% names(off)) |>
    mutate(cum_pos = win_mid + off[chr])

  gw_mean <- mean(d$fst, na.rm = TRUE)
  gw_top1 <- quantile(d$fst, 0.99, na.rm = TRUE)

  # Chromosome label tick positions
  chr_ticks <- chr_sizes |>
    mutate(tick = (offset + size / 2) / 1e6)

  # Alternating light/dark chromosome backgrounds
  bg_even <- chr_sizes |> filter(as.integer(row_number()) %% 2 == 0)

  p <- ggplot(d, aes(cum_pos / 1e6, fst)) +
    geom_rect(data = bg_even,
              aes(xmin = offset / 1e6, xmax = (offset + size) / 1e6),
              ymin = -Inf, ymax = Inf,
              fill = "grey93", inherit.aes = FALSE) +
    geom_point(size = 0.35, alpha = 0.5, colour = "grey30") +
    geom_hline(yintercept = gw_mean,
               linetype = 1, colour = "steelblue4", linewidth = 0.5) +
    geom_hline(yintercept = gw_top1,
               linetype = 2, colour = "firebrick", linewidth = 0.4) +
    scale_x_continuous(breaks = chr_ticks$tick, labels = chr_ticks$chr,
                       expand = expansion(mult = 0.01)) +
    labs(x = NULL, y = expression(F[ST]~"(Hudson, windowed)"),
         title = title,
         subtitle = sprintf(
           "%d kb windows, min %d SNPs.  Solid = genome mean (%.3f).  Dashed = top 1%% (%.3f).",
           win_size / 1000, MIN_SNPS, gw_mean, gw_top1)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.text.x        = element_text(face = "bold"))

  if (!is.null(inv_tbl)) {
    inv_d <- inv_tbl |>
      filter(chr %in% names(off)) |>
      mutate(cum_start = (start + off[chr]) / 1e6,
             cum_end   = (end   + off[chr]) / 1e6,
             label_x   = case_when(
               name == "In(3R)K" ~ cum_start,                    # left — avoids In(3R)P
               TRUE               ~ (cum_start + cum_end) / 2    # centre for all others
             ))
    p <- p +
      geom_rect(data = inv_d,
                aes(xmin = cum_start, xmax = cum_end),
                ymin = -Inf, ymax = Inf,
                fill = "steelblue", alpha = 0.18, inherit.aes = FALSE) +
      geom_text(data = inv_d,
                aes(x = label_x, label = name),
                y = Inf, vjust = 1.4, size = 2.8, colour = "steelblue4",
                inherit.aes = FALSE)
  }
  p
}

p_dm_50 <- plot_fst(dm_50,
  "Toronto vs Miami (D. mel) — windowed Fₛₜ (50 kb)",
  win_size = WIN_50, inv_tbl = inversions_mel)

p_ds_50 <- plot_fst(ds_50,
  "Tallahassee-1 vs Tallahassee-2 (D. sim) — windowed Fₛₜ (50 kb)",
  win_size = WIN_50)

p_dm_10 <- plot_fst(dm_10,
  "Toronto vs Miami (D. mel) — windowed Fₛₜ (10 kb)",
  win_size = WIN_10, inv_tbl = inversions_mel)

p_ds_10 <- plot_fst(ds_10,
  "Tallahassee-1 vs Tallahassee-2 (D. sim) — windowed Fₛₜ (10 kb)",
  win_size = WIN_10)

ggsave(file.path(out_dir, "fig_fst_Toronto_vs_Miami_50kb.png"),
       p_dm_50, width = 14, height = 4, dpi = 200)
ggsave(file.path(out_dir, "fig_fst_Tallahassee_50kb.png"),
       p_ds_50, width = 14, height = 4, dpi = 200)
ggsave(file.path(out_dir, "fig_fst_Toronto_vs_Miami_10kb.png"),
       p_dm_10, width = 14, height = 4, dpi = 200)
ggsave(file.path(out_dir, "fig_fst_Tallahassee_10kb.png"),
       p_ds_10, width = 14, height = 4, dpi = 200)

cat("\nWrote output to:", out_dir, "\n")
