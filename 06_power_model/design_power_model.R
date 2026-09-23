#!/usr/bin/env Rscript
# design_power_model.R
# -----------------------------------------------------------------------------
# Reusable model comparing detection power of two GWAS designs for a major gene:
#   (1) Truncation (artificial) selection with R replicate lines  [empirical null]
#   (2) Replicated case/control extreme-pool design               [analytic null]
#
# Everything is conditional on a gene's starting allele frequency (p0) and its
# per-allele effect size (e, in phenotypic SD). The truncation NOISE is measured
# empirically from an observed selection line; the case/control noise is analytic.
#
# Swap the lineage (or chromosome / contrast) by changing the `lines` list in the
# driver at the bottom. No parameter is hard-coded inside the functions.
# -----------------------------------------------------------------------------
suppressMessages(library(tidyverse))

# ============================ PARAMETERS =====================================
params <- list(
  # ---- case/control comparator (from the literature; no empirical data) ----
  cc_N        = 1000,             # individuals per pool
  cc_cov      = c(100, 200),      # sequencing coverage levels to show
  cc_frac     = 0.05,             # extreme fraction selected (top 5%)
  # ---- truncation selection schedule (Sarah's experiment) ----
  sel_frac    = c(0.800,0.732,0.644,0.596,0.528,0.460,0.392,0.324,
                  0.256,0.188,0.168,0.148,0.120,0.100,0.084,0.084),
  n_scored    = 500,              # flies scored per generation
  R_levels    = c(1,2,5,10),      # replicate-line counts to show
  # ---- hypothetical major gene grid ----
  p0_levels   = c(0.05,0.15,0.35),
  e_grid      = seq(0.05,1.5,0.025),
  # ---- multiple testing ----
  n_snps      = 2e6,
  alpha       = 0.05,
  freq_binw   = 0.02              # starting-freq bin width for empirical noise
)

# ============================ CORE FUNCTIONS =================================
sel_intensity <- function(f) dnorm(qnorm(1 - f)) / f

# read one selection line: returns per-SNP starting freq p0 and observed change dp
read_line <- function(refalt_file, ctrl_col, sel_col, censor = NULL) {
  d <- read.table(refalt_file, header = TRUE, sep = "", check.names = FALSE)
  af <- function(alt, ref) alt / (alt + ref)
  out <- tibble(
    chr = if ("CHROM" %in% names(d)) d$CHROM else NA,
    pos = d$POS,
    p0  = af(d[[paste0("ALT_", ctrl_col)]], d[[paste0("REF_", ctrl_col)]]),
    ps  = af(d[[paste0("ALT_", sel_col )]], d[[paste0("REF_", sel_col )]])
  ) |>
    filter(is.finite(p0), is.finite(ps), p0 > 0, p0 < 1) |>
    mutate(dp = ps - p0)
  if (!is.null(censor))
    out <- out |> filter(!(chr == censor$chrom & pos >= censor$start & pos <= censor$end))
  out
}

# empirical single-line conditional noise SD as a smooth function of p0
build_sigma_er <- function(line, binw = 0.02, min_n = 100) {
  bs <- line |>
    mutate(fb = cut(p0, breaks = seq(0, 1, binw), include.lowest = TRUE)) |>
    group_by(fb) |>
    summarise(p_mid = mean(p0), n = n(), s = sd(dp), .groups = "drop") |>
    filter(n >= min_n, s > 0)
  list(fn = approxfun(bs$p_mid, bs$s, rule = 2), bins = bs)
}

# expected truncation response of a causal allele over the whole schedule
er_signal <- function(e, p0, er_i) {
  p <- p0
  for (ig in er_i) p <- min(max(p + ig * e * p * (1 - p), 0), 1)
  p - p0
}
# case/control expected shift (single extreme-pool truncation step)
cc_signal <- function(e, p0, i_cc) i_cc * e * p0 * (1 - p0)
# case/control analytic sampling SD
cc_sigma  <- function(p, N, cov) sqrt(p * (1 - p) * (1 / N + 2 / cov))

# power at genome-wide threshold from a noncentrality parameter
make_power_fn <- function(n_snps, alpha) {
  zc <- qnorm(1 - 0.5 * (alpha / n_snps))
  function(ncp) pnorm(ncp - zc) + pnorm(-ncp - zc)
}

# ---- full power grid for one selection line ----
run_power <- function(line, params, label) {
  er_i  <- sel_intensity(round(params$sel_frac * params$n_scored) / params$n_scored)
  i_cc  <- sel_intensity(params$cc_frac)
  sig   <- build_sigma_er(line, params$freq_binw)$fn
  powr  <- make_power_fn(params$n_snps, params$alpha)
  g <- expand_grid(p0 = params$p0_levels, e = params$e_grid)
  er <- map_dfr(params$R_levels, function(R) g |>
    transmute(label, p0, e, design = paste0("Trunc R=", R),
      power = powr(sqrt(R) * map2_dbl(e, p0, \(e,p) er_signal(e,p,er_i)) / sig(p0))))
  cc <- map_dfr(params$cc_cov, function(cov) g |>
    transmute(label, p0, e, design = paste0("CC ", cov, "X"),
      power = powr(map2_dbl(e, p0, \(e,p) cc_signal(e,p,i_cc)) / cc_sigma(p0, params$cc_N, cov))))
  bind_rows(er, cc)
}

# ============================ AGGREGATION HELPERS ============================
# Pool many strata (chr x population, and in principle x direction) into one
# genome-wide empirical noise model, and propagate between-stratum uncertainty
# into power via a block bootstrap that resamples whole strata.

# load every (chromosome x population) stratum for the Dark contrast
load_all_strata <- function(dir, chroms, censor = NULL) {
  map_dfr(chroms, function(ch) {
    f <- file.path(dir, paste0("RefAlt.", ch, ".txt"))
    map_dfr(c("HOULE", "HOUSTON"), function(pop) {
      ln <- read_line(f, paste0(pop, "_L1F"), paste0(pop, "_L3F"), censor)
      ln |> mutate(pop = pop, strat = paste(ch, pop, sep = "."))
    })
  })
}

# per-stratum, per-bin sufficient statistics (so each bootstrap draw is cheap)
suff_stats <- function(dat, bw = 0.02) {
  mids <- seq(bw/2, 1 - bw/2, bw)
  dat |> mutate(bin = pmin(floor(p0/bw) + 1L, length(mids))) |>
    group_by(strat, bin) |>
    summarise(n = n(), s1 = sum(dp), s2 = sum(dp^2), .groups = "drop") |>
    mutate(p_mid = mids[bin])
}

# sigma_ER(p) from a chosen (possibly bootstrap-resampled) set of strata
sigma_from_suff <- function(suff, sel_strata, min_n = 100) {
  agg <- suff |> filter(strat %in% sel_strata) |>
    group_by(bin, p_mid) |>
    summarise(n = sum(n), s1 = sum(s1), s2 = sum(s2), .groups = "drop") |>
    filter(n >= min_n) |>
    mutate(sd = sqrt(pmax((s2 - s1^2/n)/(n - 1), 0)))
  approxfun(agg$p_mid, agg$sd, rule = 2)
}

# genome-wide power with a 90% block-bootstrap CI over strata
run_power_uncertainty <- function(dat, params, B = 400, seed = 1) {
  suff   <- suff_stats(dat, params$freq_binw)
  strata <- unique(suff$strat)
  er_i   <- sel_intensity(round(params$sel_frac * params$n_scored) / params$n_scored)
  i_cc   <- sel_intensity(params$cc_frac)
  powr   <- make_power_fn(params$n_snps, params$alpha)
  grid   <- expand_grid(p0 = params$p0_levels, e = params$e_grid) |>
    mutate(dp_er = map2_dbl(e, p0, \(e,p) er_signal(e, p, er_i)),
           dp_cc = map2_dbl(e, p0, \(e,p) cc_signal(e, p, i_cc)))
  pw <- function(sig) bind_rows(
    map_dfr(params$R_levels, \(R) grid |> transmute(p0, e, design = paste0("Trunc R=", R),
              power = powr(sqrt(R) * dp_er / sig(p0)))),
    map_dfr(params$cc_cov, \(cov) grid |> transmute(p0, e, design = paste0("CC ", cov, "X"),
              power = powr(dp_cc / cc_sigma(p0, params$cc_N, cov)))))
  central <- pw(sigma_from_suff(suff, strata)) |> rename(power_c = power)
  set.seed(seed)
  BM <- do.call(cbind, map(1:B, function(b)
    pw(sigma_from_suff(suff, sample(strata, length(strata), replace = TRUE)))$power))
  central$lo <- apply(BM, 1, quantile, 0.05)
  central$hi <- apply(BM, 1, quantile, 0.95)
  central
}

# ============================ DRIVER =========================================
if (sys.nframe() == 0) {
  dir    <- "."
  chroms <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")
  censor <- list(chrom = "chr3L", start = 8500000, end = 9000000)  # Dmel artifact region
  # Data available here: Dark-vs-Control (L3F vs L1F), both populations, 5 chromosomes.
  # The Light direction (L2F) is a second contrast to add as further strata when available.
  dat <- load_all_strata(dir, chroms, censor)
  cat(sprintf("loaded %d SNPs across %d strata\n", nrow(dat), length(unique(dat$strat))))
  central <- run_power_uncertainty(dat, params, B = 400)
  write_csv(central, "power_uncertainty.csv")
  cat("wrote power_uncertainty.csv\n")
}
