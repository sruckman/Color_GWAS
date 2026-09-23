#!/usr/bin/env Rscript
# plot_chr3L_cnv.R
#
# Three-panel figure zoomed to chr3L 8.5-9.0 Mb (Dmel coordinates):
#   Panel 1: Manhattan (-log10 p) — all 4 populations, Dsim mapped to Dmel coords
#   Panel 2: Coverage (normalized depth) — all 12 samples, Dsim mapped to Dmel coords
#   Panel 3: Delta allele frequency — all 4 populations, Dsim mapped to Dmel coords
#
# Uses raw (non-downsampled) GWAS results from Dm_qq_corrected / Ds_qq_corrected.
# One output PNG per comparison. Coverage panel is identical across all three files.
# Output: process/down_sample/bam_ds/plots_censored/chr3L_CNV_{comp}.png

library(tidyverse)
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork", repos = "https://cloud.r-project.org",
                   lib = .libPaths()[1])
}
library(patchwork)

# ============================================================
# SETTINGS
# ============================================================

zoom_chr   <- "chr3L"
zoom_start <- 8.5e6   # Dmel bp
zoom_end   <- 9.0e6   # Dmel bp

dm_gwas_dir <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/Dm_qq_corrected"
ds_gwas_dir <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/Ds_qq_corrected"
dm_cov_dir  <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/Dm_coverage"
ds_cov_dir  <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/Ds_coverage"
dm_refalt   <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/Dm_color"
ds_refalt   <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/color"
out_dir     <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comparisons <- list(
  list(name  = "Light_vs_Control",
       dm_g1 = c("HOULE_L2F",  "HOUSTON_L2F"),
       dm_g2 = c("HOULE_L1F",  "HOUSTON_L1F"),
       ds_g1 = c("Dsim_1-L2F", "Dsim_2-L2F"),
       ds_g2 = c("Dsim_1-L1F", "Dsim_2-L1F")),
  list(name  = "Dark_vs_Light",
       dm_g1 = c("HOULE_L3F",  "HOUSTON_L3F"),
       dm_g2 = c("HOULE_L2F",  "HOUSTON_L2F"),
       ds_g1 = c("Dsim_1-L3F", "Dsim_2-L3F"),
       ds_g2 = c("Dsim_1-L2F", "Dsim_2-L2F")),
  list(name  = "Dark_vs_Control",
       dm_g1 = c("HOULE_L3F",  "HOUSTON_L3F"),
       dm_g2 = c("HOULE_L1F",  "HOUSTON_L1F"),
       ds_g1 = c("Dsim_1-L3F", "Dsim_2-L3F"),
       ds_g2 = c("Dsim_1-L1F", "Dsim_2-L1F"))
)

# Population colors (matches rest of project)
pop_colors <- c(
  "D. mel-1" = "dodgerblue",
  "D. mel-2" = "darkorchid2",
  "D. sim-1" = "red4",
  "D. sim-2" = "deeppink2"
)

pop_labels <- c(
  "D. mel-1" = expression(italic("D. mel")*"-1"),
  "D. mel-2" = expression(italic("D. mel")*"-2"),
  "D. sim-1" = expression(italic("D. sim")*"-1"),
  "D. sim-2" = expression(italic("D. sim")*"-2")
)

# Coverage: 12 samples, population hue with 3 shades per timepoint (C/L/D)
cov_colors <- c(
  "HOULE_L1F"   = "#9ecae1",
  "HOULE_L2F"   = "#3182bd",
  "HOULE_L3F"   = "#08519c",
  "HOUSTON_L1F" = "#bcbddc",
  "HOUSTON_L2F" = "#756bb1",
  "HOUSTON_L3F" = "#54278f",
  "Dsim_1-L1F"  = "#fc9272",
  "Dsim_1-L2F"  = "#de2d26",
  "Dsim_1-L3F"  = "#67000d",
  "Dsim_2-L1F"  = "#fa9fb5",
  "Dsim_2-L2F"  = "#f768a1",
  "Dsim_2-L3F"  = "#ae017e"
)

cov_labels <- c(
  "HOULE_L1F"   = expression(italic("D. mel")*"-1 C"),
  "HOULE_L2F"   = expression(italic("D. mel")*"-1 L"),
  "HOULE_L3F"   = expression(italic("D. mel")*"-1 D"),
  "HOUSTON_L1F" = expression(italic("D. mel")*"-2 C"),
  "HOUSTON_L2F" = expression(italic("D. mel")*"-2 L"),
  "HOUSTON_L3F" = expression(italic("D. mel")*"-2 D"),
  "Dsim_1-L1F"  = expression(italic("D. sim")*"-1 C"),
  "Dsim_1-L2F"  = expression(italic("D. sim")*"-1 L"),
  "Dsim_1-L3F"  = expression(italic("D. sim")*"-1 D"),
  "Dsim_2-L1F"  = expression(italic("D. sim")*"-2 C"),
  "Dsim_2-L2F"  = expression(italic("D. sim")*"-2 L"),
  "Dsim_2-L3F"  = expression(italic("D. sim")*"-2 D")
)

dm_samples <- names(cov_colors)[1:6]
ds_samples <- names(cov_colors)[7:12]

# ============================================================
# COORDINATE MAPPING: Dsim chr3L -> Dmel chr3L
# Synteny segment: sim 7,851,416-8,834,090 <-> mel 8,000,000-9,000,000 (not inverted)
# ============================================================

SIM_SEG_START <- 7851416L
SIM_SEG_END   <- 8834090L
MEL_SEG_START <- 8000000L
MEL_SEG_END   <- 9000000L

sim_to_mel <- function(sim_pos) {
  MEL_SEG_START + (sim_pos - SIM_SEG_START) /
    (SIM_SEG_END - SIM_SEG_START) *
    (MEL_SEG_END - MEL_SEG_START)
}

# ============================================================
# HELPERS
# ============================================================

# Raw GWAS: Dm_qq_corrected / Ds_qq_corrected, tab-delimited, column LOG10_P_corrected
read_gwas_raw <- function(gwas_dir, pop, comp, chr) {
  f <- file.path(gwas_dir, paste0("GWAS_corrected_", pop, "_", comp, "_", chr, ".txt"))
  if (!file.exists(f)) { cat(sprintf("  [missing] %s\n", basename(f))); return(NULL) }
  read_delim(f, delim = "\t", show_col_types = FALSE) %>%
    mutate(CHROM = chr, POP = pop)
}

read_coverage <- function(cov_dir, sample, chr) {
  f <- file.path(cov_dir, paste0(sample, "_", chr, "_1kb.txt"))
  if (!file.exists(f)) { cat(sprintf("  [missing cov] %s\n", basename(f))); return(NULL) }
  df <- read.table(f, header = TRUE, sep = "\t")
  chr_mean <- mean(df$MEAN_DEPTH, na.rm = TRUE)
  if (is.na(chr_mean) || chr_mean == 0) return(NULL)
  df$NORM_DEPTH <- df$MEAN_DEPTH / chr_mean
  df$SAMPLE     <- sample
  df
}

freq_safe <- function(alt, ref, min_cov = 10) {
  total <- alt + ref
  ifelse(total >= min_cov, alt / total, NA_real_)
}

calc_delta <- function(refalt_dir, g1_samples, g2_samples, pop_label, chr) {
  f <- file.path(refalt_dir, paste0("RefAlt.", chr, ".txt"))
  if (!file.exists(f)) { cat(sprintf("  [missing refalt] %s\n", basename(f))); return(NULL) }
  df <- read.table(f, header = TRUE, sep = "", check.names = FALSE)
  p_g1 <- rowMeans(sapply(g1_samples, function(s) {
    freq_safe(as.numeric(df[[paste0("ALT_", s)]]),
              as.numeric(df[[paste0("REF_", s)]]))
  }), na.rm = TRUE)
  p_g2 <- rowMeans(sapply(g2_samples, function(s) {
    freq_safe(as.numeric(df[[paste0("ALT_", s)]]),
              as.numeric(df[[paste0("REF_", s)]]))
  }), na.rm = TRUE)
  data.frame(POS = df$POS, DELTA = abs(p_g1 - p_g2), POP = pop_label) %>%
    filter(!is.na(DELTA), is.finite(DELTA))
}

# ============================================================
# COVERAGE PANEL (same for all comparisons — build once)
# ============================================================

cat("Reading coverage...\n")

cov_dm <- map_dfr(dm_samples, ~ read_coverage(dm_cov_dir, .x, zoom_chr)) %>%
  filter(!is.na(NORM_DEPTH),
         WIN_START >= zoom_start, WIN_START <= zoom_end) %>%
  rename(MEL_POS = WIN_START)

cov_ds <- map_dfr(ds_samples, ~ read_coverage(ds_cov_dir, .x, zoom_chr)) %>%
  filter(!is.na(NORM_DEPTH),
         WIN_START >= SIM_SEG_START, WIN_START <= SIM_SEG_END) %>%
  mutate(MEL_POS = sim_to_mel(WIN_START)) %>%
  filter(MEL_POS >= zoom_start, MEL_POS <= zoom_end)

cov_all <- bind_rows(cov_dm, cov_ds) %>%
  mutate(SAMPLE = factor(SAMPLE, levels = names(cov_colors)))

cat(sprintf("  Dmel cov rows: %d | Dsim cov rows: %d\n", nrow(cov_dm), nrow(cov_ds)))

p_cov <- ggplot(cov_all, aes(x = MEL_POS / 1e6, y = NORM_DEPTH,
                              color = SAMPLE, group = SAMPLE)) +
  geom_line(linewidth = 0.4, alpha = 0.8) +
  scale_color_manual(values = cov_colors, labels = cov_labels, name = NULL) +
  scale_x_continuous(limits = c(zoom_start / 1e6, zoom_end / 1e6)) +
  coord_cartesian(ylim = c(0, 8)) +
  labs(x = NULL, y = "Normalized\ncoverage") +
  theme_bw(base_size = 11) +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.4, "cm"),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

# ============================================================
# MAIN LOOP — one output file per comparison
# ============================================================

for (comp in comparisons) {
  cat(sprintf("\n=== %s ===\n", comp$name))

  comp_label <- gsub("_", " ", comp$name)

  # ---- Manhattan ----
  cat("  Building Manhattan...\n")

  gwas_dm <- map_dfr(c("HOULE", "HOUSTON"), function(pop) {
    read_gwas_raw(dm_gwas_dir, pop, comp$name, zoom_chr)
  }) %>%
    mutate(POP     = recode(POP, "HOULE" = "D. mel-1", "HOUSTON" = "D. mel-2"),
           MEL_POS = POS) %>%
    filter(!is.na(LOG10_P_corrected), is.finite(LOG10_P_corrected),
           LOG10_P_corrected >= 0,
           MEL_POS >= zoom_start, MEL_POS <= zoom_end)

  gwas_ds <- map_dfr(c("Dsim_1", "Dsim_2"), function(pop) {
    read_gwas_raw(ds_gwas_dir, pop, comp$name, zoom_chr)
  }) %>%
    filter(!is.na(LOG10_P_corrected), is.finite(LOG10_P_corrected),
           LOG10_P_corrected >= 0,
           POS >= SIM_SEG_START, POS <= SIM_SEG_END) %>%
    mutate(MEL_POS = sim_to_mel(POS),
           POP     = recode(POP, "Dsim_1" = "D. sim-1", "Dsim_2" = "D. sim-2")) %>%
    filter(MEL_POS >= zoom_start, MEL_POS <= zoom_end)

  gwas_all <- bind_rows(gwas_dm, gwas_ds) %>%
    mutate(POP = factor(POP, levels = names(pop_colors)))

  cat(sprintf("    Dmel SNPs: %d | Dsim SNPs: %d\n", nrow(gwas_dm), nrow(gwas_ds)))

  p_man <- ggplot(gwas_all,
                  aes(x = MEL_POS / 1e6, y = LOG10_P_corrected,
                      color = POP, alpha = LOG10_P_corrected)) +
    geom_point(size = 0.3) +
    scale_color_manual(values = pop_colors, labels = pop_labels, name = NULL) +
    scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
    scale_x_continuous(limits = c(zoom_start / 1e6, zoom_end / 1e6)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = NULL,
         y = expression(-log[10](italic(p))),
         title = paste0("chr3L CNV region — ", comp_label)) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "right",
          legend.text      = element_text(size = 9),
          panel.grid.minor = element_blank(),
          axis.text.x      = element_blank(),
          axis.ticks.x     = element_blank())

  # ---- Delta frequency ----
  cat("  Building delta freq...\n")

  delta_dm <- bind_rows(
    calc_delta(dm_refalt, comp$dm_g1[grepl("HOULE",   comp$dm_g1)],
               comp$dm_g2[grepl("HOULE",   comp$dm_g2)], "D. mel-1", zoom_chr),
    calc_delta(dm_refalt, comp$dm_g1[grepl("HOUSTON", comp$dm_g1)],
               comp$dm_g2[grepl("HOUSTON", comp$dm_g2)], "D. mel-2", zoom_chr)
  ) %>%
    mutate(MEL_POS = POS) %>%
    filter(MEL_POS >= zoom_start, MEL_POS <= zoom_end)

  delta_ds <- bind_rows(
    calc_delta(ds_refalt,
               comp$ds_g1[grepl("Dsim_1", comp$ds_g1)],
               comp$ds_g2[grepl("Dsim_1", comp$ds_g2)], "D. sim-1", zoom_chr),
    calc_delta(ds_refalt,
               comp$ds_g1[grepl("Dsim_2", comp$ds_g1)],
               comp$ds_g2[grepl("Dsim_2", comp$ds_g2)], "D. sim-2", zoom_chr)
  ) %>%
    filter(POS >= SIM_SEG_START, POS <= SIM_SEG_END) %>%
    mutate(MEL_POS = sim_to_mel(POS)) %>%
    filter(MEL_POS >= zoom_start, MEL_POS <= zoom_end)

  delta_all <- bind_rows(delta_dm, delta_ds) %>%
    mutate(POP = factor(POP, levels = names(pop_colors)))

  cat(sprintf("    Dmel delta SNPs: %d | Dsim delta SNPs: %d\n",
              nrow(delta_dm), nrow(delta_ds)))

  p_delta <- ggplot(delta_all,
                    aes(x = MEL_POS / 1e6, y = DELTA,
                        color = POP, alpha = DELTA)) +
    geom_point(size = 0.3) +
    scale_color_manual(values = pop_colors, labels = pop_labels, name = NULL) +
    scale_alpha_continuous(range = c(0.2, 0.7), guide = "none") +
    scale_x_continuous(limits = c(zoom_start / 1e6, zoom_end / 1e6)) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "chr3L position (Mb, D. mel coords)",
         y = expression("|" ~ Delta ~ "allele freq|")) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "right",
          legend.text      = element_text(size = 9),
          panel.grid.minor = element_blank())

  # ---- Combine and save ----
  p_combined <- p_man / p_cov / p_delta +
    plot_layout(heights = c(1.2, 1, 1), guides = "keep")

  outfile <- file.path(out_dir, paste0("chr3L_CNV_", comp$name, ".png"))
  ggsave(outfile, p_combined, width = 9, height = 10, dpi = 150)
  cat(sprintf("  Saved: %s\n", basename(outfile)))

  rm(gwas_dm, gwas_ds, gwas_all, delta_dm, delta_ds, delta_all,
     p_man, p_delta, p_combined)
  gc()
}

cat("\n--- Done ---\n")
