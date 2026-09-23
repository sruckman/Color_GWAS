#!/usr/bin/env Rscript
# make_figures_Ds.R
#
# D. simulans wrapper around Tony Long's make_figures.R.
# Sources design_power_model.R, overrides load_all_strata() for Dsim paths
# and sample names, then runs all five figures for both contrasts.
#
# Figures written to:
#   plots_censored/design_power_Dsim_{comp}/fig_*.png
#
# Run from /dfs7/adl/sruckman/XQTL/XQTL2 via run_make_figures_Ds.sh

suppressMessages(library(tidyverse))
suppressMessages(library(ggtext))
source("scripts/design_power_model.R")   # core functions + params
params$n_snps <- 3.5e6                  # Dsim-specific: ~3.5M SNPs tested per contrast

# ---- Override load_all_strata for Dsim ----------------------------------------
# Dsim sample names use hyphens: Dsim_1-L1F, Dsim_1-L2F, Dsim_1-L3F etc.
load_all_strata <- function(dir, chroms, censor = NULL, sel_tp = "L3F") {
  pops <- list(
    list(label = "Tallahassee-1", code = "Dsim_1"),
    list(label = "Tallahassee-2", code = "Dsim_2")
  )
  map_dfr(chroms, function(ch) {
    f <- file.path(dir, paste0("RefAlt.", ch, ".txt"))
    if (!file.exists(f)) { cat(sprintf("  [missing] %s\n", basename(f))); return(NULL) }
    map_dfr(pops, function(p) {
      ln <- read_line(f,
                      paste0(p$code, "-L1F"),
                      paste0(p$code, "-", sel_tp),
                      NULL)   # pass NULL; apply multi-window censor below
      if (!is.null(censor))
        for (w in censor)
          ln <- ln |> filter(!(chr == w$chrom & pos >= w$start & pos <= w$end))
      ln |> mutate(pop = p$label, strat = paste(ch, p$label, sep = "."))
    })
  })
}

# ============================================================
# SETTINGS
# ============================================================
ds_base  <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"
out_base <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored"

CHROMS <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
CENSOR <- list(
  list(chrom = "chr3L", start = 8342753L, end = 8834090L),
  list(chrom = "chrX",  start = 7800000L, end = 8400000L)   # Confirmed by coverage scan
)
B_BOOT <- 400
M_NULL <- 200000

# shared palette / theme (Tony's make_figures.R)
lev <- c("Trunc R=1","Trunc R=2","Trunc R=5","Trunc R=10","Trunc R=20","CC 100X","CC 200X")
pal <- c("Trunc R=1"="#7a0177","Trunc R=2"="#c51b8a","Trunc R=5"="#f768a1",
         "Trunc R=10"="#fb9a99","Trunc R=20"="#d9a0a0","CC 100X"="#08519c","CC 200X"="#6baed6")
theme_pub <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title    = element_text(size = 12, face = "plain"),
        plot.subtitle = element_markdown(size = 8.4))

# ============================================================
# RUN FOR EACH CONTRAST
# ============================================================
contrasts <- list(
  list(name = "Light_vs_Control", sel_tp = "L2F"),
  list(name = "Dark_vs_Control",  sel_tp = "L3F")
)

for (ctr in contrasts) {
  comp   <- ctr$name
  SEL    <- ctr$sel_tp
  DIR    <- file.path(ds_base, comp, "Ds_color_sep")
  OUTDIR <- file.path(out_base, paste0("design_power_Dsim_", comp))
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  fp   <- function(f) file.path(OUTDIR, f)

  SCOPE <- sprintf(
    "*D. simulans* (Tallahassee-1 + Tallahassee-2) — %s — %d chromosomes",
    gsub("_", " ", comp), length(CHROMS))
  sub2 <- function(detail) paste0(SCOPE, "\n", detail)

  cat(sprintf("\n=== Dsim %s ===\n", comp))
  set.seed(1)

  # Load genome-wide data
  dat <- load_all_strata(DIR, CHROMS, CENSOR, sel_tp = SEL)
  cat(sprintf("  %d SNPs across %d strata\n", nrow(dat), length(unique(dat$strat))))

  # Shared model objects
  er_i <- sel_intensity(round(params$sel_frac * params$n_scored) / params$n_scored)
  i_cc <- sel_intensity(params$cc_frac)
  powr <- make_power_fn(params$n_snps, params$alpha)
  zc   <- qnorm(1 - 0.5 * (params$alpha / params$n_snps))

  suff_gw  <- suff_stats(dat, params$freq_binw)
  strata   <- unique(suff_gw$strat)
  sig_gw   <- sigma_from_suff(suff_gw, strata)
  sig_bins <- suff_gw |>
    group_by(bin, p_mid) |>
    summarise(n = sum(n), s1 = sum(s1), s2 = sum(s2), .groups = "drop") |>
    filter(n >= 100) |>
    mutate(s = sqrt(pmax((s2 - s1^2/n) / (n - 1), 0))) |>
    arrange(p_mid)

  # ------------------------------------------------------------------
  # FIG A  fig_conditional_sd.png
  # ------------------------------------------------------------------
  message("  A  fig_conditional_sd ...")
  condsd <- bind_rows(
    map_dfr(c(1,2,5,10,20), \(R) tibble(
      design = paste0("R=", R), kind = "Truncation selection (empirical)",
      p0 = sig_bins$p_mid, sd = sig_bins$s / sqrt(R))),
    map_dfr(c(100,200), \(cov) tibble(
      design = paste0(cov, "X"), kind = "Case/control (modeled)",
      p0 = sig_bins$p_mid, sd = cc_sigma(sig_bins$p_mid, params$cc_N, cov)))
  )
  condsd$design <- factor(condsd$design,
    levels = c("R=1","R=2","R=5","R=10","R=20","100X","200X"))
  csd_pal <- c("R=1"="#7a0177","R=2"="#c51b8a","R=5"="#f768a1",
               "R=10"="#fbb4b9","R=20"="#fde0dd","100X"="#08519c","200X"="#6baed6")
  fA <- ggplot(condsd, aes(p0, sd, colour = design, linetype = kind)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.1) +
    scale_colour_manual(values = csd_pal, name = NULL) +
    scale_linetype_manual(
      values = c("Truncation selection (empirical)" = "solid",
                 "Case/control (modeled)" = "22"), name = NULL) +
    labs(x = "Starting (control) allele frequency",
         y = "SD of allele-frequency change",
         title = "Background noise depends on starting frequency, not just design",
         subtitle = sub2("Truncation (solid, empirical) vs modeled case/control (dashed).")) +
    theme_pub + theme(legend.position = "right")
  ggsave(fp("fig_conditional_sd.png"), fA, width = 9, height = 5.5, dpi = 200)

  # ------------------------------------------------------------------
  # FIG B  fig_signal_vs_snr.png
  # ------------------------------------------------------------------
  message("  B  fig_signal_vs_snr ...")
  snr <- expand_grid(p0 = c(0.05, 0.15, 0.25, 0.35), e = 1.0) |>
    mutate(dp_trunc = map2_dbl(e, p0, \(e,p) er_signal(e, p, er_i)),
           dp_cc    = map2_dbl(e, p0, \(e,p) cc_signal(e, p, i_cc)),
           z_trunc  = dp_trunc / sig_gw(p0),
           z_cc     = dp_cc / cc_sigma(p0, params$cc_N, 100))
  plotd <- snr |>
    select(p0, dp_trunc, dp_cc, z_trunc, z_cc) |>
    pivot_longer(-p0) |>
    mutate(
      panel  = ifelse(grepl("^dp", name),
                      "Absolute shift |Δp|  (what intuition tracks)",
                      "Signal-to-noise z  (what detection needs)"),
      design = ifelse(grepl("trunc", name), "Truncation (1 line)", "Case/control 100X"))
  fB <- ggplot(plotd, aes(p0, value, colour = design)) +
    geom_line(linewidth = 1) + geom_point(size = 2.4) +
    facet_wrap(~panel, scales = "free_y") +
    scale_colour_manual(
      values = c("Truncation (1 line)" = "#7a0177",
                 "Case/control 100X"   = "#08519c"), name = NULL) +
    geom_hline(
      data = tibble(panel = "Signal-to-noise z  (what detection needs)", y = zc),
      aes(yintercept = y), linetype = 3, colour = "grey55") +
    labs(x = "Starting allele frequency p0", y = NULL,
         title = "Truncation selection moves a major allele more, yet detects it less",
         subtitle = sub2("Large-effect allele (e = 1 SD). Left: absolute shift. Right: signal-to-noise (dotted = genome-wide z threshold).")) +
    theme_pub
  ggsave(fp("fig_signal_vs_snr.png"), fB, width = 9.5, height = 4.8, dpi = 200)

  # ------------------------------------------------------------------
  # FIG C  fig_null_qq.png
  # ------------------------------------------------------------------
  message("  C  fig_null_qq (Monte Carlo) ...")
  binw  <- params$freq_binw
  resid <- dat |>
    mutate(fb = pmin(floor(p0/binw) + 1L, round(1/binw))) |>
    group_by(fb) |> mutate(r = dp - mean(dp)) |> ungroup()
  idx    <- sample(nrow(dat), M_NULL, replace = TRUE)
  p_null <- dat$p0[idx]; s_null <- sig_gw(p_null)
  resid_by_bin <- split(resid$r, resid$fb)
  draw_R <- function(R) {
    fb <- pmin(floor(p_null/binw) + 1L, round(1/binw))
    m  <- numeric(M_NULL)
    for (b in unique(fb)) {
      sel  <- which(fb == b)
      pool <- resid_by_bin[[as.character(b)]]
      if (is.null(pool)) pool <- resid$r
      m[sel] <- rowMeans(matrix(sample(pool, length(sel)*R, replace = TRUE), ncol = R))
    }
    m
  }
  L10 <- function(z) -log10(pmax(2 * pnorm(-abs(z)), .Machine$double.xmin))
  allL10 <- list()
  for (R in c(1,2,5,10,20))
    allL10[[paste0("Trunc R=", R)]] <- L10(draw_R(R) / (s_null / sqrt(R)))
  for (cov in c(100,200)) {
    dpn <- rnorm(M_NULL, 0, cc_sigma(p_null, params$cc_N, cov))
    allL10[[paste0("CC ", cov, "X")]] <- L10(dpn / cc_sigma(p_null, params$cc_N, cov))
  }
  qq <- imap_dfr(allL10, function(v, nm) {
    v <- sort(v); n <- length(v)
    keep <- unique(round(seq(1, n, length.out = 3000)))
    tibble(design = nm, expected = -log10((n - keep + 0.5)/n), observed = v[keep])
  })
  qq$design <- factor(qq$design, levels = lev)
  fC <- ggplot(filter(qq, expected <= 3.3), aes(expected, observed, colour = design)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = pal, name = NULL) +
    coord_cartesian(xlim = c(0,3.3), ylim = c(0,3.3)) +
    labs(x = "Expected -log10(p)  (uniform null)",
         y = "Observed -log10(p)",
         title = "Null -log10(p) is well-calibrated; truncation R=1 is conservative",
         subtitle = sub2("On grey line = calibrated. R=1 is conservative: one line's change is bounded.")) +
    theme_pub + theme(legend.position = "right")
  ggsave(fp("fig_null_qq.png"), fC, width = 9, height = 5.5, dpi = 200)

  # ------------------------------------------------------------------
  # FIG D  fig_power.png
  # ------------------------------------------------------------------
  message("  D  fig_power (bootstrap) ...")
  central <- run_power_uncertainty(dat, params, B = B_BOOT, seed = 1)
  write_csv(central, fp("power_uncertainty.csv"))
  cp <- central |>
    mutate(design = factor(design, levels = lev),
           pf     = factor(paste0("p0 = ", p0)))
  fD <- ggplot(cp, aes(e, power_c, colour = design, fill = design)) +
    geom_hline(yintercept = 0.8, linetype = 3, colour = "grey60") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~pf) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values   = pal, name = NULL) +
    scale_y_continuous(limits = c(0,1)) +
    labs(x = "Per-allele effect size e (phenotypic SD)",
         y = "Detection power (genome-wide)",
         title = "Power to detect a major gene, integrated over genome-wide noise uncertainty",
         subtitle = sub2("Bands = 90% block-bootstrap CI over chr x pop strata. Case/control has no lineage input.")) +
    theme_pub + theme(legend.position = "right", plot.title = element_text(size = 11.5))
  ggsave(fp("fig_power.png"), fD, width = 10, height = 4.8, dpi = 200)

  # ------------------------------------------------------------------
  # FIG E  fig_design_space.png
  # ------------------------------------------------------------------
  message("  E  fig_design_space ...")
  G0 <- expand_grid(p0 = seq(0.02, 0.5, 0.002), e = seq(0.05, 1.5, 0.008)) |>
    mutate(dp_er = map2_dbl(e, p0, \(e,p) er_signal(e, p, er_i)),
           dp_cc = map2_dbl(e, p0, \(e,p) cc_signal(e, p, i_cc)),
           pw_cc = powr(dp_cc / cc_sigma(p0, params$cc_N, 200)))
  Rlev <- c(2,5,10); thr <- 0.8
  G <- map_dfr(Rlev, \(R) G0 |>
    mutate(R = R,
           pw_tr = powr(sqrt(R) * dp_er / sig_gw(p0)),
           region = case_when(
             pw_tr >= thr & pw_cc >= thr ~ "Both detect",
             pw_tr >= thr & pw_cc <  thr ~ "Truncation only",
             pw_tr <  thr & pw_cc >= thr ~ "Case/control only",
             TRUE                        ~ "Neither detects")))
  G$region <- factor(G$region,
    levels = c("Neither detects","Case/control only","Truncation only","Both detect"))
  G$panel  <- factor(paste0("Truncation R = ", G$R),
                     levels = paste0("Truncation R = ", Rlev))
  Bb  <- 200
  bnd <- map_dfr(Rlev, function(R) {
    ev <- G0 |> distinct(p0, e, dp_er); fr <- numeric(nrow(ev))
    for (b in seq_len(Bb)) {
      sg <- sigma_from_suff(suff_gw, sample(strata, length(strata), replace = TRUE))
      fr <- fr + (powr(sqrt(R) * ev$dp_er / sg(ev$p0)) >= thr)
    }
    ev |> mutate(R = R, detect_frac = fr / Bb)
  })
  bnd$panel <- factor(paste0("Truncation R = ", bnd$R), levels = levels(G$panel))
  cols <- c("Neither detects"   = "grey86",
            "Case/control only" = "#6baed6",
            "Truncation only"   = "#c51b8a",
            "Both detect"       = "#7a3f9e")
  fE <- ggplot(G, aes(p0, e, fill = region)) +
    geom_raster() +
    geom_contour(data = bnd, aes(p0, e, z = detect_frac),
                 breaks = c(0.05, 0.95), colour = "grey20",
                 linewidth = 0.35, linetype = "22", inherit.aes = FALSE) +
    facet_wrap(~panel) +
    scale_fill_manual(values = cols, name = NULL) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Starting allele frequency p0",
         y = "Effect size e (phenotypic SD)",
         title = "Where each design wins, and how replication expands the truncation advantage",
         subtitle = sub2("Truncation R lines vs replicated case/control (N=1000, 200X); 80% power. Dashed = 5-95% bootstrap band.")) +
    theme_bw(base_size = 12) +
    theme(panel.grid     = element_blank(), legend.position = "right",
          plot.title     = element_text(size = 10.5),
          plot.subtitle  = element_markdown(size = 7.9),
          panel.spacing  = unit(0.9, "lines"))
  ggsave(fp("fig_design_space.png"), fE, width = 12, height = 4.4, dpi = 200)

  cat(sprintf("  Done. Figures in %s\n", OUTDIR))
  rm(dat, fA, fB, fC, fD, fE, central, suff_gw, sig_gw, sig_bins,
     G0, G, bnd, condsd, snr, plotd, qq, allL10, resid); gc()
}

cat("\n--- Done ---\n")
