#!/usr/bin/env Rscript
# plot_chorion_evidence.R
#
# 3-panel figure showing the chorion amplification artifact at chr3L and chrX.
# Run separately for each species (Dmel, Dsim) and each contrast.
#
# Panel layout — facet_grid(panel ~ chr):
#   Row A  -log10(p) from Fisher GWAS, every 5th SNP
#   Row B  Normalized 1 kb coverage (genome-wide median = 1; dashed line)
#   Row C  |Δp| between selected and control pool, 5 kb window means
#
# Columns: chr3L (left) and chrX (right).
# Red shading marks the chorion cluster windows in each column.
#
# Output: plots_censored/chorion_evidence/chorion_evidence_{sp}_{comp}.png

suppressMessages(library(tidyverse))

# ============================================================
# SETTINGS
# ============================================================
dm_base    <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"
ds_base    <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"
dm_cov_dir <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/Dm_coverage"
ds_cov_dir <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/Ds_coverage"
out_dir    <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/chorion_evidence"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_chroms   <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
focus_chroms <- c("chr3L", "chrX")

# Chorion windows — NATIVE coordinates per species.
# D. mel: 3L 66D (8.5-9.0 Mb) and X 7F (8.4-8.6 Mb; Cp36 at 8,479,319)
# D. sim: 3L 66D (8.34-8.83 Mb) and X 7F (8.70-8.78 Mb PLACEHOLDER — see note below)
# NOTE: Dsim X coordinates are approximate. Synteny from tan (Dsim 8,770,795) predicts
# the real 7F cluster at ~8.03 Mb. The current window (8.70-8.78 Mb) likely covers tan.
# Rerun check_chorion_clusters_Ds.R with chrX bracket 7.0-8.5 Mb to get real coordinates.
chorion_mel <- tribble(
  ~chr,    ~start,   ~end,
  "chr3L", 8500000L, 9000000L,
  "chrX",  8400000L, 8600000L
)
chorion_sim <- tribble(
  ~chr,    ~start,   ~end,
  "chr3L", 8342753L, 8834090L,
  "chrX",  7800000L, 8400000L   # Confirmed by coverage scan. Peak at 8,208,001 in
                                 # 5/6 lines (up to 68x); secondary peak 7,860,001.
                                 # Old 8.7 Mb window was tan, not chorion.
)

pop_colors_dm <- c(HOULE   = "dodgerblue",  HOUSTON = "darkorchid2")
pop_colors_ds <- c(Dsim_1  = "red4",        Dsim_2  = "deeppink2")
pop_labels_dm <- c(HOULE   = "Toronto",     HOUSTON = "Miami")
pop_labels_ds <- c(Dsim_1  = "Tallahassee-1", Dsim_2 = "Tallahassee-2")

PANEL_GWAS <- "-log10(p)"
PANEL_COV  <- "Normalized\ncoverage"
PANEL_DP   <- "|Δp| (5 kb mean)"
panel_levels <- c(PANEL_GWAS, PANEL_COV, PANEL_DP)

GWAS_THIN <- 1L     # keep every SNP (zoomed region is small enough)
DP_WIN    <- 5000L  # |Δp| binning window (bp)

# Zoom windows — tight view around each chorion cluster (per species).
# The facet x-axis will be set by the data within these limits.
# Dmel chrX extends to 9.8 Mb to include tan at 9.22 Mb.
# Dsim chrX is very tight (8.4-9.2 Mb) because the 7F cluster is only 82 kb wide;
#   tan at 8.770,795 falls inside/adjacent to the chorion window (8.702-8.784 Mb).
zoom_mel <- tribble(
  ~chr,    ~lo,     ~hi,
  "chr3L", 7.5e6,   10.5e6,   # 66D cluster 8.5-9.0 Mb
  "chrX",  7.5e6,    9.8e6    # 7F cluster 8.4-8.6 Mb; tan at 9.22 Mb visible
)
zoom_sim <- tribble(
  ~chr,    ~lo,     ~hi,
  "chr3L", 7.5e6,   10.0e6,   # 66D cluster 8.34-8.83 Mb
  "chrX",  7.0e6,    9.5e6    # broad view: predicted real 7F cluster ~8.03 Mb + tan at 8.77 Mb
)

# tan gene positions for chrX annotation (vertical line + label on chrX panel only)
# Dmel: tan at 9,217,655 — clearly outside the 7F chorion window (8.4-8.6 Mb)
# Dsim: tan at 8,770,795 — shown to confirm it is NOT in the chorion cluster;
#   the broad zoom (7.0-9.5 Mb) will show both the predicted cluster (~8.03 Mb)
#   and tan (~8.77 Mb), letting you see they are NOT the same feature.
tan_mel <- tibble(chr = "chrX", pos = 9217655L, label = "tan")
tan_sim <- tibble(chr = "chrX", pos = 8770795L, label = "tan")

comparisons <- list(
  list(name = "Light_vs_Control", ctrl_tp = "L1F", sel_tp = "L2F"),
  list(name = "Dark_vs_Control",  ctrl_tp = "L1F", sel_tp = "L3F")
)

# ============================================================
# FUNCTIONS
# ============================================================

# Load GWAS -log10(p) for chr3L and chrX, thinned for plotting.
load_gwas_panel <- function(gwas_dir, pops, comp,
                             chrs = focus_chroms, thin = GWAS_THIN) {
  map_dfr(chrs, function(ch) {
    map_dfr(pops, function(pop) {
      f <- file.path(gwas_dir,
                     paste0("GWAS_results_", pop, "_", comp, "_", ch, ".txt"))
      if (!file.exists(f)) { cat(sprintf("    [missing] %s\n", basename(f))); return(NULL) }
      read_delim(f, delim = " ", show_col_types = FALSE) |>
        filter(!is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0) |>
        mutate(chr = ch, pop = pop)
    })
  }) |>
    group_by(chr, pop) |>
    arrange(POS, .by_group = TRUE) |>
    slice(seq(1L, n(), thin)) |>
    ungroup() |>
    transmute(chr, pos = POS, pop, sample = pop, value = LOG10_P,
              panel = factor(PANEL_GWAS, levels = panel_levels))
}

# Load 1 kb coverage files; normalize each sample to its genome-wide median.
# All chromosomes are loaded for normalization, then filtered to chrs_plot.
load_coverage_panel <- function(cov_dir, samples, chrs_plot = focus_chroms) {
  map_dfr(samples, function(samp) {
    all_d <- map_dfr(all_chroms, function(ch) {
      f <- file.path(cov_dir, paste0(samp, "_", ch, "_1kb.txt"))
      if (!file.exists(f)) return(NULL)
      d <- read.table(f, header = TRUE)
      tibble(chr = ch, pos = d$WIN_START, depth = d$MEAN_DEPTH)
    })
    if (nrow(all_d) == 0) {
      cat(sprintf("    [no coverage] %s\n", samp)); return(NULL)
    }
    gw_med <- median(all_d$depth[all_d$depth > 0], na.rm = TRUE)
    all_d |>
      filter(chr %in% chrs_plot, !is.na(depth)) |>
      mutate(norm = depth / gw_med, sample = samp,
             pop  = sub("[-_](L[123]F)$", "", samp))
  }) |>
    transmute(chr, pos, pop, sample, value = norm,
              panel = factor(PANEL_COV, levels = panel_levels))
}

# Compute |Δp| from RefAlt files, binned to DP_WIN-bp windows.
# sep: "_" for Dmel (HOULE_L1F), "-" for Dsim (Dsim_1-L1F).
load_dp_panel <- function(refalt_dir, pops, ctrl_tp, sel_tp,
                           sep = "_", chrs = focus_chroms, win = DP_WIN) {
  map_dfr(chrs, function(ch) {
    f <- file.path(refalt_dir, paste0("RefAlt.", ch, ".txt"))
    if (!file.exists(f)) { cat(sprintf("    [missing] RefAlt.%s.txt\n", ch)); return(NULL) }
    d <- read.table(f, header = TRUE, sep = "", check.names = FALSE) |> as_tibble()

    map_dfr(pops, function(pop) {
      ctrl_id <- paste0(pop, sep, ctrl_tp)
      sel_id  <- paste0(pop, sep, sel_tp)
      rc  <- d[[paste0("REF_", ctrl_id)]];  ac  <- d[[paste0("ALT_", ctrl_id)]]
      rs  <- d[[paste0("REF_", sel_id)]];   as_ <- d[[paste0("ALT_", sel_id)]]
      if (is.null(rc) || is.null(rs)) {
        warning("Columns not found for ", ctrl_id, " / ", sel_id, " in ", ch)
        return(NULL)
      }
      p_ctrl <- ac  / (rc  + ac)
      p_sel  <- as_ / (rs  + as_)
      tibble(chr = ch, pos = d$POS, pop = pop,
             abs_dp = abs(p_sel - p_ctrl)) |>
        filter(is.finite(abs_dp))
    })
  }) |>
    mutate(win_start = floor((pos - 1L) / win) * win) |>
    group_by(chr, pop, win_start) |>
    summarise(pos = as.integer(win_start + win / 2L),
              value = mean(abs_dp), .groups = "drop") |>
    transmute(chr, pos, pop, sample = pop, value,
              panel = factor(PANEL_DP, levels = panel_levels))
}

# Build the 3-panel figure from a combined long tibble.
make_figure <- function(all_d, chorion_tbl, pop_cols, pop_labs,
                        species_label, comp_label, gene_marks = NULL) {
  all_d <- all_d |>
    mutate(chr   = factor(chr,   levels = focus_chroms),
           panel = factor(panel, levels = panel_levels))

  rects <- chorion_tbl |>
    mutate(chr = factor(chr, levels = focus_chroms))

  # Horizontal reference line at 1 in coverage panel only
  cov_ref <- all_d |>
    filter(panel == PANEL_COV) |>
    distinct(panel, chr) |>
    mutate(yint = 1)

  ggplot() +
    # Chorion shading — appears in all panel rows for each chr column
    geom_rect(
      data = rects,
      aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf),
      fill = "firebrick", alpha = 0.12, inherit.aes = FALSE
    ) +
    # Panel A: GWAS -log10(p) as points
    geom_point(
      data        = filter(all_d, panel == PANEL_GWAS),
      aes(x = pos / 1e6, y = value, colour = pop),
      size = 0.25, alpha = 0.5
    ) +
    # Panel B: normalized coverage — one line per sample
    geom_line(
      data        = filter(all_d, panel == PANEL_COV),
      aes(x = pos / 1e6, y = value, colour = pop, group = sample),
      linewidth = 0.25, alpha = 0.65
    ) +
    geom_hline(
      data = cov_ref, aes(yintercept = yint),
      linetype = 2, colour = "grey55", linewidth = 0.4
    ) +
    # Panel C: |Δp| 5 kb means — one line per population
    geom_line(
      data        = filter(all_d, panel == PANEL_DP),
      aes(x = pos / 1e6, y = value, colour = pop, group = pop),
      linewidth = 0.8
    ) +
    # Gene marker lines (e.g. tan on chrX) — all panels, labelled in top panel only
    { if (!is.null(gene_marks)) list(
        geom_vline(
          data = gene_marks %>% mutate(chr = factor(chr, levels = focus_chroms)),
          aes(xintercept = pos / 1e6),
          colour = "black", linetype = "dotted", linewidth = 0.5,
          inherit.aes = FALSE
        ),
        geom_label(
          data = gene_marks %>%
            mutate(chr = factor(chr, levels = focus_chroms),
                   panel = factor(PANEL_GWAS, levels = panel_levels)),
          aes(x = pos / 1e6, y = Inf, label = label),
          vjust = 1.3, hjust = 0.5, size = 3,
          colour = "black", fill = "white", label.size = 0,
          label.padding = unit(0.1, "lines"), inherit.aes = FALSE
        )
      ) else list() } +
    scale_colour_manual(values = pop_cols, labels = pop_labs, name = NULL) +
    facet_grid(panel ~ chr, scales = "free") +
    labs(
      x        = "Position (Mb)",
      y        = NULL,
      title    = paste0(species_label, " — chorion cluster amplification artifact"),
      subtitle = comp_label
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text.y     = element_text(angle = 0),
      strip.text.x     = element_text(face = "bold"),
      legend.position  = "top",
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(0.55, "lines"),
      plot.subtitle    = element_text(size = 9, colour = "grey35")
    )
}

# ============================================================
# PRE-LOAD COVERAGE
# Coverage files are not comparison-specific, so load once.
# ============================================================
dm_samples <- c("HOULE_L1F",  "HOULE_L2F",  "HOULE_L3F",
                 "HOUSTON_L1F","HOUSTON_L2F","HOUSTON_L3F")
ds_samples <- c("Dsim_1-L1F","Dsim_1-L2F","Dsim_1-L3F",
                 "Dsim_2-L1F","Dsim_2-L2F","Dsim_2-L3F")

dm_pops <- c("HOULE", "HOUSTON")
ds_pops <- c("Dsim_1", "Dsim_2")

cat("Loading D. mel coverage (all chromosomes for genome-wide normalization)...\n")
dm_cov <- load_coverage_panel(dm_cov_dir, dm_samples)
cat(sprintf("  %d coverage rows\n", nrow(dm_cov)))

cat("Loading D. sim coverage...\n")
ds_cov <- load_coverage_panel(ds_cov_dir, ds_samples)
cat(sprintf("  %d coverage rows\n", nrow(ds_cov)))

# ============================================================
# MAIN LOOP — one figure per species per contrast
# ============================================================
for (ctr in comparisons) {
  comp   <- ctr$name
  ctrl_t <- ctr$ctrl_tp
  sel_t  <- ctr$sel_tp
  clabel <- gsub("_", " ", comp)

  cat(sprintf("\n=== %s ===\n", comp))

  dm_dir <- file.path(dm_base, comp, "Dm_color_sep")
  ds_dir <- file.path(ds_base, comp, "Ds_color_sep")

  # ---- D. melanogaster ----
  cat("  D. mel: GWAS and |Δp|...\n")
  dm_gwas <- load_gwas_panel(dm_dir, dm_pops, comp)
  dm_dp   <- load_dp_panel(dm_dir, dm_pops, ctrl_t, sel_t, sep = "_")
  cat(sprintf("    GWAS rows: %d   |Δp| rows: %d\n", nrow(dm_gwas), nrow(dm_dp)))

  # Apply zoom: filter each dataset to the tight window around each cluster.
  # facet_grid with scales = "free" will then set per-column x limits from the data.
  apply_zoom <- function(d, zoom_tbl) {
    d |>
      inner_join(zoom_tbl, by = "chr") |>
      filter(pos >= lo, pos <= hi) |>
      select(-lo, -hi)
  }

  dm_all <- bind_rows(dm_gwas, dm_cov, dm_dp) |> apply_zoom(zoom_mel)
  if (nrow(dm_all) > 0) {
    p <- make_figure(dm_all, chorion_mel, pop_colors_dm, pop_labels_dm,
                     "D. melanogaster", clabel, gene_marks = tan_mel)
    outf <- file.path(out_dir, sprintf("chorion_evidence_Dmel_%s.png", comp))
    ggsave(outf, p, width = 10, height = 9, dpi = 200)
    cat(sprintf("  Saved %s\n", basename(outf)))
    rm(p)
  }
  rm(dm_gwas, dm_dp, dm_all); gc()

  # ---- D. simulans ----
  cat("  D. sim: GWAS and |Δp|...\n")
  ds_gwas <- load_gwas_panel(ds_dir, ds_pops, comp)
  ds_dp   <- load_dp_panel(ds_dir, ds_pops, ctrl_t, sel_t, sep = "-")
  cat(sprintf("    GWAS rows: %d   |Δp| rows: %d\n", nrow(ds_gwas), nrow(ds_dp)))

  ds_all <- bind_rows(ds_gwas, ds_cov, ds_dp) |> apply_zoom(zoom_sim)
  if (nrow(ds_all) > 0) {
    p <- make_figure(ds_all, chorion_sim, pop_colors_ds, pop_labels_ds,
                     "D. simulans", clabel, gene_marks = tan_sim)
    outf <- file.path(out_dir, sprintf("chorion_evidence_Dsim_%s.png", comp))
    ggsave(outf, p, width = 10, height = 9, dpi = 200)
    cat(sprintf("  Saved %s\n", basename(outf)))
    rm(p)
  }
  rm(ds_gwas, ds_dp, ds_all); gc()
}

cat("\nDone. Figures in:", out_dir, "\n")
