#!/usr/bin/env Rscript
# run_design_power_Ds.R
#
# Runs Tony Long's design power model on D. simulans RefAlt data.
# Identical logic to run_design_power_Dm.R; differs in paths, sample names,
# censor region, and species label.
#
# Output: process/down_sample/bam_ds/plots_censored/
#   fig_power_Dsim_{comp}.png
#   fig_design_space_Dsim_{comp}.png
#   power_Dsim_{comp}.csv

suppressMessages(library(tidyverse))

# ============================================================
# PARAMETERS
# ============================================================
params <- list(
  cc_N      = 1000,
  cc_cov    = c(100, 200),
  cc_frac   = 0.05,
  sel_frac  = c(0.800, 0.732, 0.644, 0.596, 0.528, 0.460, 0.392, 0.324,
                0.256, 0.188, 0.168, 0.148, 0.120, 0.100, 0.084, 0.084),
  n_scored  = 500,
  R_levels  = c(1, 2, 5, 10),
  p0_levels = c(0.05, 0.15, 0.35),
  e_grid    = seq(0.05, 1.5, 0.025),
  n_snps    = 3.5e6,
  alpha     = 0.05,
  freq_binw = 0.02
)

ds_base <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"
out_dir <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chroms <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
censor <- list(
  list(chrom = "chr3L", start = 8342753L, end = 8834090L),
  list(chrom = "chrX",  start = 7800000L, end = 8400000L)   # Confirmed by coverage scan (peak at 8,208,001)
)

spp_label <- "D. simulans"  # used in plot titles (italicised via bquote)

# ============================================================
# CORE FUNCTIONS (from design_power_model.R, Long lab 2025)
# ============================================================
sel_intensity <- function(f) dnorm(qnorm(1 - f)) / f

read_line <- function(refalt_file, ctrl_col, sel_col, censor = NULL) {
  d <- read.table(refalt_file, header = TRUE, sep = "", check.names = FALSE)
  af <- function(alt, ref) alt / (alt + ref)
  out <- tibble(
    chr = if ("CHROM" %in% names(d)) d$CHROM else NA_character_,
    pos = d$POS,
    p0  = af(d[[paste0("ALT_", ctrl_col)]], d[[paste0("REF_", ctrl_col)]]),
    ps  = af(d[[paste0("ALT_", sel_col )]], d[[paste0("REF_", sel_col )]])
  ) |>
    filter(is.finite(p0), is.finite(ps), p0 > 0, p0 < 1) |>
    mutate(dp = ps - p0)
  if (!is.null(censor))
    for (w in censor)
      out <- out |> filter(!(chr == w$chrom & pos >= w$start & pos <= w$end))
  out
}

er_signal <- function(e, p0, er_i) {
  p <- p0
  for (ig in er_i) p <- min(max(p + ig * e * p * (1 - p), 0), 1)
  p - p0
}

cc_signal <- function(e, p0, i_cc) i_cc * e * p0 * (1 - p0)

cc_sigma <- function(p, N, cov) sqrt(p * (1 - p) * (1 / N + 2 / cov))

make_power_fn <- function(n_snps, alpha) {
  zc <- qnorm(1 - 0.5 * (alpha / n_snps))
  function(ncp) pnorm(ncp - zc) + pnorm(-ncp - zc)
}

suff_stats <- function(dat, bw = 0.02) {
  mids <- seq(bw / 2, 1 - bw / 2, bw)
  dat |>
    mutate(bin = pmin(floor(p0 / bw) + 1L, length(mids))) |>
    group_by(strat, bin) |>
    summarise(n = n(), s1 = sum(dp), s2 = sum(dp^2), .groups = "drop") |>
    mutate(p_mid = mids[bin])
}

sigma_from_suff <- function(suff, sel_strata, min_n = 100) {
  agg <- suff |>
    filter(strat %in% sel_strata) |>
    group_by(bin, p_mid) |>
    summarise(n = sum(n), s1 = sum(s1), s2 = sum(s2), .groups = "drop") |>
    filter(n >= min_n) |>
    mutate(sd = sqrt(pmax((s2 - s1^2 / n) / (n - 1), 0)))
  approxfun(agg$p_mid, agg$sd, rule = 2)
}

run_power_uncertainty <- function(dat, params, B = 400, seed = 1) {
  suff   <- suff_stats(dat, params$freq_binw)
  strata <- unique(suff$strat)
  er_i   <- sel_intensity(round(params$sel_frac * params$n_scored) / params$n_scored)
  i_cc   <- sel_intensity(params$cc_frac)
  powr   <- make_power_fn(params$n_snps, params$alpha)

  grid <- expand_grid(p0 = params$p0_levels, e = params$e_grid) |>
    mutate(dp_er = map2_dbl(e, p0, \(e, p) er_signal(e, p, er_i)),
           dp_cc = map2_dbl(e, p0, \(e, p) cc_signal(e, p, i_cc)))

  pw <- function(sig) bind_rows(
    map_dfr(params$R_levels, \(R) grid |>
              transmute(p0, e, design = paste0("E&R  R=", R),
                        power = powr(sqrt(R) * dp_er / sig(p0)))),
    map_dfr(params$cc_cov,   \(cov) grid |>
              transmute(p0, e, design = paste0("CC  ", cov, "X"),
                        power = powr(dp_cc / cc_sigma(p0, params$cc_N, cov))))
  )

  central <- pw(sigma_from_suff(suff, strata)) |> rename(power_c = power)
  set.seed(seed)
  BM <- do.call(cbind, map(seq_len(B), function(b)
    pw(sigma_from_suff(suff,
                       sample(strata, length(strata), replace = TRUE)))$power))
  central$lo <- apply(BM, 1, quantile, 0.05)
  central$hi <- apply(BM, 1, quantile, 0.95)
  central
}

run_design_space <- function(dat, params,
                              p0_grid = seq(0.02, 0.50, 0.02),
                              e_dense = seq(0.05, 1.5, 0.025),
                              B = 200, seed = 2) {
  suff   <- suff_stats(dat, params$freq_binw)
  strata <- unique(suff$strat)
  er_i   <- sel_intensity(round(params$sel_frac * params$n_scored) / params$n_scored)
  i_cc   <- sel_intensity(params$cc_frac)
  powr   <- make_power_fn(params$n_snps, params$alpha)

  grid <- expand_grid(p0 = p0_grid, e = e_dense) |>
    mutate(dp_er = map2_dbl(e, p0, \(e, p) er_signal(e, p, er_i)),
           dp_cc = map2_dbl(e, p0, \(e, p) cc_signal(e, p, i_cc)))

  pow_cc <- powr(grid$dp_cc / cc_sigma(grid$p0, params$cc_N, max(params$cc_cov)))

  er_pow_central <- map_dfr(params$R_levels, function(R) {
    grid |> transmute(p0, e, R,
                      pow_er = powr(sqrt(R) * dp_er / sigma_from_suff(suff, strata)(p0)),
                      pow_cc = pow_cc)
  })

  set.seed(seed)
  boot_bounds <- map_dfr(seq_len(B), function(b) {
    sig_b <- sigma_from_suff(suff, sample(strata, length(strata), replace = TRUE))
    map_dfr(params$R_levels, function(R) {
      pw <- powr(sqrt(R) * grid$dp_er / sig_b(grid$p0))
      tibble(p0 = grid$p0, e = grid$e, R = R, pow = pw) |>
        group_by(p0, R) |>
        summarise(e80 = min(e[pow >= 0.8], default = NA_real_), .groups = "drop") |>
        mutate(boot = b)
    })
  })

  list(central = er_pow_central, boot_bounds = boot_bounds)
}

# ============================================================
# DATA LOADING
# Dsim sample names use hyphens: Dsim_1-L1F, Dsim_1-L2F, Dsim_1-L3F
# ============================================================
load_all_strata_Ds <- function(comp_name, chroms, censor = NULL) {
  sel_tp   <- if (grepl("Light", comp_name)) "L2F" else "L3F"
  comp_dir <- file.path(ds_base, comp_name, "Ds_color_sep")

  pops <- list(
    list(label = "D. sim-1", code = "Dsim_1"),
    list(label = "D. sim-2", code = "Dsim_2")
  )

  map_dfr(chroms, function(ch) {
    f <- file.path(comp_dir, paste0("RefAlt.", ch, ".txt"))
    if (!file.exists(f)) {
      cat(sprintf("  [missing] %s\n", basename(f)))
      return(NULL)
    }
    map_dfr(pops, function(p) {
      # Dsim sample names use hyphens: e.g. Dsim_1-L1F
      ln <- read_line(f,
                      paste0(p$code, "-L1F"),
                      paste0(p$code, "-", sel_tp),
                      censor)
      ln |> mutate(pop = p$label, strat = paste(ch, p$label, sep = "."))
    })
  })
}

# ============================================================
# PLOTTING
# ============================================================
design_colours <- c(
  "E&R  R=1"   = "#bdd7ee",
  "E&R  R=2"   = "#5b9bd5",
  "E&R  R=5"   = "#2e6da4",
  "E&R  R=10"  = "#1a3f6b",
  "CC  100X"   = "#f4b183",
  "CC  200X"   = "#c55a11"
)

p0_labs <- c("0.05" = "p₀ = 0.05",
             "0.15" = "p₀ = 0.15",
             "0.35" = "p₀ = 0.35")

plot_power <- function(pw, comp_label, spp) {
  pw |>
    mutate(p0 = factor(p0)) |>
    ggplot(aes(e, power_c, colour = design, fill = design)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0.8, linetype = 2, colour = "grey40") +
    facet_wrap(~p0, labeller = as_labeller(p0_labs)) +
    scale_colour_manual(values = design_colours, name = "Design") +
    scale_fill_manual(values   = design_colours, name = "Design") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(1)) +
    labs(
      x        = "Per-allele effect size (phenotypic SD)",
      y        = "Detection power",
      title    = bquote(italic(.(spp)) ~ "—" ~ .(gsub("_", " ", comp_label))),
      subtitle = "Dashed line = 80% power.  Bands = 90% bootstrap CI (block-resampled strata)."
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "bottom",
          panel.grid.minor = element_blank())
}

plot_design_space <- function(ds, comp_label, spp, R_show = c(2, 5, 10)) {
  R_labels <- paste0("E&R  R = ", R_show)

  cells <- ds$central |>
    filter(R %in% R_show) |>
    mutate(
      wins = case_when(
        pow_er >= 0.8 & pow_cc >= 0.8 ~ "Both",
        pow_er >= 0.8                  ~ "E&R only",
        pow_cc >= 0.8                  ~ "CC only",
        TRUE                           ~ "Neither"
      ),
      wins    = factor(wins, levels = c("Neither", "CC only", "E&R only", "Both")),
      R_label = factor(paste0("E&R  R = ", R), levels = R_labels)
    )

  boot_sum <- ds$boot_bounds |>
    filter(R %in% R_show) |>
    group_by(p0, R) |>
    summarise(e80_lo = quantile(e80, 0.05, na.rm = TRUE),
              e80_hi = quantile(e80, 0.95, na.rm = TRUE),
              .groups = "drop") |>
    mutate(R_label = factor(paste0("E&R  R = ", R), levels = R_labels))

  dp <- diff(sort(unique(cells$p0)))[1]
  de <- diff(sort(unique(cells$e)))[1]

  ggplot(cells, aes(p0, e)) +
    geom_tile(aes(fill = wins), width = dp, height = de) +
    geom_ribbon(data = boot_sum,
                aes(x = p0, ymin = e80_lo, ymax = e80_hi),
                inherit.aes = FALSE, fill = "white", alpha = 0.40) +
    facet_wrap(~R_label) +
    scale_fill_manual(
      values = c("Neither"  = "grey88",
                 "CC only"  = "#f4b183",
                 "E&R only" = "#5b9bd5",
                 "Both"     = "#70ad47"),
      name = "≥80% power"
    ) +
    labs(
      x        = "Starting allele frequency (p₀)",
      y        = "Per-allele effect size (SD)",
      title    = bquote(italic(.(spp)) ~ "—" ~ .(gsub("_", " ", comp_label))),
      subtitle = paste0("CC: N=", params$cc_N, ", ", max(params$cc_cov),
                        "X.  White band = 90% bootstrap CI on E&R 80%-power boundary.")
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid      = element_blank())
}

# ============================================================
# RUN
# ============================================================
for (comp in c("Light_vs_Control", "Dark_vs_Control")) {
  cat(sprintf("\n=== Dsim %s ===\n", comp))

  dat <- load_all_strata_Ds(comp, chroms, censor)
  cat(sprintf("  %d SNPs across %d strata\n", nrow(dat), length(unique(dat$strat))))

  ## Power curves with bootstrap CI
  pw <- run_power_uncertainty(dat, params, B = 400)
  write_csv(pw, file.path(out_dir, paste0("power_Dsim_", comp, ".csv")))

  ggsave(file.path(out_dir, paste0("fig_power_Dsim_", comp, ".png")),
         plot_power(pw, comp, spp_label),
         width = 9, height = 5, dpi = 150)
  cat("  Saved: fig_power\n")

  ## Design-space map
  ds <- run_design_space(dat, params, B = 200)
  ggsave(file.path(out_dir, paste0("fig_design_space_Dsim_", comp, ".png")),
         plot_design_space(ds, comp, spp_label),
         width = 9, height = 5, dpi = 150)
  cat("  Saved: fig_design_space\n")

  rm(dat, pw, ds); gc()
}

cat("\n--- Done ---\n")
