#!/usr/bin/env Rscript
# simulate_null_qq.R
#
# Compares observed QQ plots against two simulated scenarios:
#   1. Drift only  — two populations diverging from the same base by drift alone
#   2. Variable polygenic selection — control line drifts; selected line has
#      drift + per-SNP selection coefficients drawn from a Gamma(shape=2, mean=s_mean).
#
# Why variable (not uniform) effects?
#   Uniform selection gives every SNP the same Δp = p*(1-p)*i*s, which is
#   bounded by p*(1-p). This compresses test-statistic variance and produces a
#   concave-down QQ curve. Drawing effects from a Gamma(shape=2) gives a few
#   SNPs large effects and many SNPs small effects, creating spread in the
#   test-statistic distribution and a linear QQ slope that matches observations.
#   With truncating selection on a highly polygenic trait, LD drags the whole
#   genome — every SNP shifts, not just QTNs — so assigning nonzero effects
#   to all SNPs is biologically appropriate.
#
# Tune s_mean to match the observed slope:
#   larger s_mean -> steeper green envelope
#
# Output: process/down_sample/bam_ds/plots_censored/QQ_null_simulation_{comp}.png

library(tidyverse)

# ============================================================
# SETTINGS
# ============================================================

N_start    <- 500    # starting census population size
N_end      <- 42     # breeders at end (8.4% of 500)
n_gen      <- 16     # generations of selection
coverage   <- 150    # simulated pool-seq depth per pool

# Exact selection proportions per generation (from experimental design figure)
sel_schedule <- c(0.800, 0.732, 0.644, 0.596, 0.528, 0.460, 0.392,
                  0.324, 0.256, 0.188, 0.168, 0.148, 0.120, 0.100,
                  0.084, 0.084)

n_snps     <- 1000000 # background SNPs per simulation (match scale of ~2M observed)
s_mean     <- 0.08    # mean of Gamma(shape=2) distribution for per-block effects;
                      # larger = STEEPER/higher green envelope. Raised from 0.04
                      # because the green line was too shallow vs the observed
                      # lines. Push toward 0.10-0.12 if still shallow; back off if
                      # the top of the green curve flattens (SNPs saturating at 0.999).
                      # Raising ld_decay/block_size steepens the upper tail too.
n_reps     <- 50      # simulation replicates
set.seed(42)

# Optional command-line override of s_mean, so several values can run in parallel:
#   Rscript simulate_null_qq.R 0.10
# Output filenames are tagged with the value, so runs do not overwrite each other.
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1 && !is.na(suppressWarnings(as.numeric(.args[1])))) {
  s_mean <- as.numeric(.args[1])
}
s_tag <- sprintf("s%03d", round(s_mean * 1000))   # 0.08 -> s080, 0.10 -> s100
cat(sprintf("Using s_mean = %.4f (output tag %s)\n", s_mean, s_tag))

# ---- LD block structure ----
# The real data has large LD blocks and inversions, so the response is smeared
# across linked SNPs rather than localized; matching the observed QQ required
# many SNPs AND substantial LD. We model LD explicitly: SNPs are grouped into
# contiguous blocks, and within a block each SNP's allele-frequency trajectory
# is a blend of (a) a shared block trajectory and (b) its own independent
# trajectory, with the weight on the shared trajectory decaying exponentially
# with distance from the block start: w = exp(-pos_in_block / ld_decay).
# This gives pairwise LD ~ exp(-d/ld_decay) within a block and ~0 between blocks
# (max LD radius = block_size). The block carries the selected variant; a SNP
# feels selection in proportion to its LD with the block.
# Three knobs to tune to the observed slope: ld_decay and block_size (more/longer
# LD = heavier tail) and s_mean. Increase ld_decay/block_size for "more LD".
block_size <- 100   # SNPs per LD block (max LD radius; ~100 kb at this SNP density)
ld_decay   <- 40    # exponential LD decay constant, in SNPs

# Paths
dm_base  <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"
ds_base  <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"
out_dir  <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

comparisons <- c("Light_vs_Control", "Dark_vs_Light", "Dark_vs_Control")
chroms      <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

censor_mel <- list(chrom = "chr3L", start = 8500000, end = 9000000)
censor_sim <- list(chrom = "chr3L", start = 8342753, end = 8834090)

pop_colors <- c(
  "D. mel-1" = "dodgerblue",
  "D. mel-2" = "darkorchid2",
  "D. sim-1" = "red4",
  "D. sim-2" = "deeppink2"
)

# ============================================================
# SCHEDULES
# ============================================================

# Ne decreases linearly as truncation intensifies
Ne_schedule <- round(seq(N_start, N_end, length.out = n_gen))

# Selection intensity i(p) = phi(x_p) / p  (Falconer & Mackay 1996)
sel_intensity <- function(p) { x <- qnorm(1 - p); dnorm(x) / p }
i_schedule <- sel_intensity(sel_schedule)

cat("Ne schedule:     ", paste(Ne_schedule,               collapse = " -> "), "\n")
cat("Sel proportion:  ", paste(round(sel_schedule, 3),    collapse = " -> "), "\n")
cat("Sel intensity:   ", paste(round(i_schedule,   3),    collapse = " -> "), "\n\n")

# ---- LD block indexing + blend helper ----
pos_in_block <- (seq_len(n_snps) - 1L) %% block_size        # 0 .. block_size-1
block_id     <- (seq_len(n_snps) - 1L) %/% block_size + 1L   # 1 .. n_blocks
n_blocks     <- max(block_id)
w_shared     <- exp(-pos_in_block / ld_decay)               # weight on shared block trajectory
cat(sprintf("LD blocks: %d blocks of %d SNPs; decay = %d SNPs\n\n",
            n_blocks, block_size, ld_decay))

# Build one population's per-SNP frequencies by blending a block-level trajectory
# (expanded to its member SNPs) with an independent per-SNP trajectory.
blend_line <- function(block_traj, indep_traj) {
  p <- w_shared * block_traj[block_id] + (1 - w_shared) * indep_traj
  pmax(0.001, pmin(0.999, p))
}

# ============================================================
# SIMULATION FUNCTIONS
# ============================================================

# Vectorised chi-square approximation to Fisher's exact test
chisq_pvals <- function(a, b, c, d) {
  n     <- a + b + c + d
  valid <- (a + b) > 0 & (c + d) > 0 & (a + c) > 0 & (b + d) > 0 & n > 0
  chi2  <- rep(NA_real_, length(a))
  chi2[valid] <- (a[valid] * d[valid] - b[valid] * c[valid])^2 * n[valid] /
    ((a[valid] + b[valid]) * (c[valid] + d[valid]) *
     (a[valid] + c[valid]) * (b[valid] + d[valid]))
  pchisq(chi2, df = 1, lower.tail = FALSE)
}

# Wright-Fisher drift only (control line)
sim_wf_drift <- function(p_init, Ne_schedule) {
  p <- p_init
  for (g in seq_along(Ne_schedule)) {
    p <- rbinom(length(p), 2L * Ne_schedule[g], p) / (2L * Ne_schedule[g])
    p <- pmax(0.001, pmin(0.999, p))
  }
  p
}

# Wright-Fisher with variable polygenic selection on all SNPs (selected line).
# Each SNP receives its own effect s_i ~ Gamma(shape=2, mean=s_mean), drawn once per replicate.
# Additive nudge per generation: Δp = p*(1-p)*i(t)*s_i
# Variable effects create spread in the chi-square distribution → linear QQ slope.
# With s_mean ~ 0.02-0.06 the max Δp per generation for average SNPs is small,
# but a tail of large-effect SNPs drives the top of the QQ.
sim_wf_var_sel <- function(p_init, Ne_schedule, i_schedule, s_effects) {
  p <- p_init
  for (g in seq_along(Ne_schedule)) {
    p <- p + p * (1 - p) * i_schedule[g] * s_effects
    p <- pmax(0.001, pmin(0.999, p))          # clip BEFORE rbinom (large s_i can push p > 1)
    p <- rbinom(length(p), 2L * Ne_schedule[g], p) / (2L * Ne_schedule[g])
    p <- pmax(0.001, pmin(0.999, p))
  }
  p
}

# Pool-seq p-values from two simulated populations.
# USE_FISHER = TRUE uses Fisher's exact test, matching the empirical pipeline
#   (run_fisher_*_bam_ds.R). This is the intended default for consistency.
#   fisher.test is not vectorised, so at n_snps = 1e6 and n_reps = 50 this is
#   slow (tens of millions of tests). Do NOT reduce n_snps to compensate; the
#   large count is needed to match the data. Instead give the job a long wall
#   time (run_null_qq.sh is set to 3 days; standard partition allows up to 14).
#   USE_FISHER = FALSE (chi-square, near-identical at coverage = 150) is only a
#   fallback for quick tests.
USE_FISHER <- TRUE

pool_pvals <- function(p_ctrl, p_sel, cov) {
  rc <- rbinom(length(p_ctrl), cov, p_ctrl)
  rs <- rbinom(length(p_sel),  cov, p_sel)
  if (USE_FISHER) {
    pv <- vapply(seq_along(rc), function(i) {
      fisher.test(matrix(c(rc[i], cov - rc[i],
                           rs[i], cov - rs[i]),
                         nrow = 2, byrow = TRUE))$p.value
    }, numeric(1))
  } else {
    pv <- chisq_pvals(rc, cov - rc, rs, cov - rs)
  }
  pv[!is.na(pv) & pv > 0]
}

# ============================================================
# RUN SIMULATIONS
# ============================================================

# Starting allele frequencies (5-95% as in Vlachos & Kofler 2019), two layers:
#   p0_block — one per LD block; carries the shared (and selected) trajectory
#   p0_snp   — one per SNP;       independent background
# Drawn once; reps differ by fresh drift draws and fresh per-block effects.
p0_block <- runif(n_blocks, 0.05, 0.95)
p0_snp   <- runif(n_snps,   0.05, 0.95)

cat(sprintf("Simulating drift-only (%d reps)...\n", n_reps))
drift_list <- vector("list", n_reps)
for (i in seq_len(n_reps)) {
  if (i %% 25 == 0) cat(sprintf("  rep %d / %d\n", i, n_reps))
  lineA <- blend_line(sim_wf_drift(p0_block, Ne_schedule),
                      sim_wf_drift(p0_snp,   Ne_schedule))
  lineB <- blend_line(sim_wf_drift(p0_block, Ne_schedule),
                      sim_wf_drift(p0_snp,   Ne_schedule))
  pv <- pool_pvals(lineA, lineB, coverage)
  drift_list[[i]] <- sort(-log10(pv), decreasing = TRUE)
}

cat(sprintf("\nSimulating variable polygenic selection (%d reps, s_mean=%.4f)...\n",
            n_reps, s_mean))
sel_list <- vector("list", n_reps)
for (i in seq_len(n_reps)) {
  if (i %% 10 == 0) cat(sprintf("  rep %d / %d\n", i, n_reps))
  # One effect per BLOCK, drawn fresh each rep (Gamma shape=2, mean=s_mean).
  # The block carries the selected variant; member SNPs feel it in proportion
  # to their LD with the block (w_shared). The independent per-SNP layer is
  # neutral background (drift only), so a SNP weakly linked to the block shows
  # mostly drift, mimicking the smearing of signal across large LD blocks.
  s_block <- rgamma(n_blocks, shape = 2, rate = 2 / s_mean)
  ctrl <- blend_line(sim_wf_drift(p0_block, Ne_schedule),
                     sim_wf_drift(p0_snp,   Ne_schedule))
  sel  <- blend_line(sim_wf_var_sel(p0_block, Ne_schedule, i_schedule, s_block),
                     sim_wf_drift(p0_snp,   Ne_schedule))
  pv <- pool_pvals(ctrl, sel, coverage)
  sel_list[[i]] <- sort(-log10(pv), decreasing = TRUE)
}

# ============================================================
# BUILD ENVELOPES (5th-95th percentile across reps)
# ============================================================

n_grid   <- 2000
# Upper bound = max expected -log10(p) for n_snps SNPs.
# E.g., n_snps = 1e6 -> x_max_sim = 6.  Setting the grid beyond this would
# cause approx() to extrapolate (flat line = plateau artifact), so we stop
# exactly at the supported range and use rule=1 to return NA outside it.
x_max_sim <- -log10(1 / n_snps)
x_shared  <- seq(0, x_max_sim, length.out = n_grid)  # evenly spaced in -log10(p) space

# Interpolate each rep's sorted -log10(p) values onto the shared x-grid,
# then take cross-rep quantiles. This avoids the straight-line artifact that
# appears when using a probability-space grid (which gives very few points
# in the tail where x > 3).
# rule=1: return NA outside the supported range so the envelope ends cleanly
# rather than being extrapolated flat.
make_envelope <- function(sim_list, label) {
  mat <- vapply(sim_list, function(x) {
    obs <- sort(x, decreasing = TRUE)
    exp <- sort(-log10(ppoints(length(obs))), decreasing = TRUE)
    approx(exp, obs, xout = x_shared, rule = 1)$y
  }, numeric(n_grid))
  tibble(
    expected  = x_shared,
    med       = apply(mat, 1, median,   na.rm = TRUE),
    lo        = apply(mat, 1, quantile, 0.05, na.rm = TRUE),
    hi        = apply(mat, 1, quantile, 0.95, na.rm = TRUE),
    condition = label
  )
}

cat("\nBuilding envelopes...\n")
env_drift <- make_envelope(drift_list, "Drift only")
env_sel   <- make_envelope(sel_list,   "Polygenic selection")

# ============================================================
# READ OBSERVED GWAS P-VALUES
# ============================================================

read_obs_pvals <- function(base_dir, subdir, pop, comp, censor) {
  map_dfr(chroms, function(chr) {
    f <- file.path(base_dir, comp, subdir,
                   paste0("GWAS_results_", pop, "_", comp, "_", chr, ".txt"))
    if (!file.exists(f)) {
      cat(sprintf("    [missing] %s\n", basename(f)))
      return(NULL)
    }
    read_delim(f, delim = " ", show_col_types = FALSE) %>%
      mutate(CHROM = chr) %>%
      filter(!is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0) %>%
      filter(!(CHROM == censor$chrom &
               POS   >= censor$start &
               POS   <= censor$end)) %>%
      select(LOG10_P)
  }) %>% pull(LOG10_P)
}

obs_specs <- list(
  list(base = dm_base, subdir = "Dm_color_sep",
       pops = c("HOULE", "HOUSTON"), names = c("D. mel-1", "D. mel-2"),
       censor = censor_mel),
  list(base = ds_base, subdir = "Ds_color_sep",
       pops = c("Dsim_1", "Dsim_2"), names = c("D. sim-1", "D. sim-2"),
       censor = censor_sim)
)

# ============================================================
# MAKE AND SAVE QQ PLOT PER COMPARISON
# ============================================================

caption_txt <- paste0(
  "Grey = drift only | Green = variable polygenic selection  ",
  "(5th-95th pct, ", n_reps, " sims each)\n",
  sprintf("N %d->%d, %d gen, s_mean=%.4f [Gamma sh=2], LD block=%d SNPs decay=%d, cov=%dx, Fisher",
          N_start, N_end, n_gen, s_mean, block_size, ld_decay, coverage)
)

for (comp in comparisons) {
  cat(sprintf("\nBuilding QQ plot: %s\n", comp))

  obs_lines <- map_dfr(obs_specs, function(sp) {
    map_dfr(seq_along(sp$pops), function(i) {
      lp <- read_obs_pvals(sp$base, sp$subdir, sp$pops[i], comp, sp$censor)
      if (length(lp) == 0) return(NULL)
      obs  <- sort(lp, decreasing = TRUE)
      exp  <- sort(-log10(ppoints(length(obs))), decreasing = TRUE)
      high <- exp > 1
      low  <- !high
      keep <- sort(c(which(high), which(low)[seq(1, sum(low), by = 50)]))
      cat(sprintf("  %s: %d SNPs (slope = %.2f)\n",
                  sp$names[i], length(obs),
                  sum(obs * exp) / sum(exp^2)))
      tibble(expected = exp[keep], observed = obs[keep], pop = sp$names[i])
    })
  })

  if (nrow(obs_lines) == 0) {
    cat("  [skipped] no observed data\n")
    next
  }

  p <- ggplot() +
    geom_ribbon(data = env_sel,
                aes(x = expected, ymin = lo, ymax = hi),
                fill = "darkgreen", alpha = 0.20) +
    geom_line(data = env_sel,
              aes(x = expected, y = med),
              color = "darkgreen", linewidth = 0.7, linetype = "dashed") +
    geom_ribbon(data = env_drift,
                aes(x = expected, ymin = lo, ymax = hi),
                fill = "grey60", alpha = 0.25) +
    geom_line(data = env_drift,
              aes(x = expected, y = med),
              color = "grey40", linewidth = 0.7, linetype = "dashed") +
    geom_line(data = obs_lines,
              aes(x = expected, y = observed, color = pop),
              linewidth = 0.6, alpha = 0.9) +
    scale_color_manual(values = pop_colors, name = NULL) +
    coord_cartesian(xlim = c(0, 7), ylim = c(0, NA)) +
    labs(
      x       = expression("Expected" ~ -log[10](italic(p))),
      y       = expression("Observed" ~ -log[10](italic(p))),
      title   = paste0("QQ — ", gsub("_", " ", comp), " (censored)"),
      caption = caption_txt
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 7, color = "grey40")
    )

  outfile <- file.path(out_dir, paste0("QQ_null_simulation_", comp, "_", s_tag, ".png"))
  ggsave(outfile, p, width = 7, height = 7, dpi = 150)
  cat(sprintf("Saved: %s\n", outfile))
}

cat("\n--- Done ---\n")
