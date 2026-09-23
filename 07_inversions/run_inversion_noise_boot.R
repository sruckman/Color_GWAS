#!/usr/bin/env Rscript
# run_inversion_noise_boot.R
#
# Block-bootstrap CIs and p-values for the inversion-noise comparison.
# Results section: "Inversions do not inflate the background noise."
#
# Miami (HOUSTON) segregates for In(2L)t, In(3L)P, In(3R)K.
# Toronto (HOULE) carries none.
#
# TEST 1 — BETWEEN POPULATIONS (autosomes only):
#   Δσ = σ(Miami) − σ(Toronto), pooled across autosomes.
#   Bootstrap unit: chromosome arm. Resample {chr2L, chr2R, chr3L, chr3R}
#   with replacement (same design as the power-model CIs in
#   run_design_power_Dm.R). Both populations contribute from the same
#   resampled arm set, preserving the pairing.
#
# TEST 2 — WITHIN POPULATION, BY ARM:
#   σ(2L)−σ(2R), σ(3L)−σ(2R), σ(3R)−σ(2R), per population.
#   Bootstrap unit: 500 kb genomic block within each arm, resampled
#   independently for each arm. Toronto is the negative control.
#
# Bootstrap p-value: 2 × min(P(Δ ≤ 0), P(Δ ≥ 0)).
# 90% CI from 400 replicates.
#
# Outputs:
#   plots_censored/inversion_noise/sigma_between_pop_boot.csv
#   plots_censored/inversion_noise/sigma_by_arm_boot.csv

suppressMessages(library(tidyverse))

# ============================================================
# PARAMETERS
# ============================================================
dm_base    <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"
out_dir    <- file.path(dm_base, "plots_censored", "inversion_noise")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FREQ_BW    <- 0.02        # frequency bin width (must match inversion_noise script)
MIN_N      <- 100         # minimum SNPs per bin to include
B          <- 400         # bootstrap replicates
CONF       <- 0.90        # CI level
P0_LEVELS  <- c(0.15, 0.35)
BLOCK_SIZE <- 500000L     # 500 kb genomic blocks for within-arm bootstrap
SEED       <- 42

auto_chroms   <- c("chr2L", "chr2R", "chr3L", "chr3R")  # X excluded throughout
inverted_arms <- c("chr2L", "chr3L", "chr3R")
ref_arm       <- "chr2R"

comps <- c("Dark_vs_Control", "Light_vs_Control")

# Chorion windows to exclude (autosomes only; chrX excluded entirely)
chorion_windows <- list(
  list(chr = "chr3L", start = 8500000L, end = 9000000L)
)

# ============================================================
# FUNCTIONS
# ============================================================
read_line <- function(refalt_file, ctrl_col, sel_col) {
  d  <- read.table(refalt_file, header = TRUE, sep = "", check.names = FALSE)
  af <- function(alt, ref) alt / (alt + ref)
  tibble(
    chr = d$CHROM,
    pos = as.integer(d$POS),
    p0  = af(d[[paste0("ALT_", ctrl_col)]], d[[paste0("REF_", ctrl_col)]]),
    ps  = af(d[[paste0("ALT_", sel_col )]], d[[paste0("REF_", sel_col )]])
  ) |>
    filter(is.finite(p0), is.finite(ps), p0 > 0, p0 < 1) |>
    mutate(dp = ps - p0)
}

apply_censor <- function(dat) {
  for (w in chorion_windows)
    dat <- filter(dat, !(chr == w$chr & pos >= w$start & pos <= w$end))
  dat
}

# Aggregate sufficient statistics for σ estimation.
# group_vars: character vector of grouping columns beyond (bin, p_mid).
make_suff <- function(dat, group_vars, bw = FREQ_BW) {
  mids <- seq(bw / 2, 1 - bw / 2, bw)
  dat |>
    mutate(bin   = pmin(floor(p0 / bw) + 1L, length(mids)),
           p_mid = mids[bin]) |>
    group_by(across(all_of(c(group_vars, "bin", "p_mid")))) |>
    summarise(n  = n(),
              s1 = sum(dp),
              s2 = sum(dp^2),
              .groups = "drop")
}

# σ(p0) from a tibble of sufficient stats (any grouping above bin level
# has already been selected). Returns NA if not enough data.
sigma_at_p0 <- function(suff_rows, p0_val, min_n = MIN_N) {
  agg <- suff_rows |>
    group_by(bin, p_mid) |>
    summarise(n = sum(n), s1 = sum(s1), s2 = sum(s2), .groups = "drop") |>
    filter(n >= min_n) |>
    mutate(sd = sqrt(pmax((s2 - s1^2 / n) / (n - 1), 0)))
  if (nrow(agg) == 0) return(NA_real_)
  approx(agg$p_mid, agg$sd, xout = p0_val, rule = 2)$y
}

load_dat <- function(comp) {
  sel_tp   <- if (grepl("Light", comp)) "L2F" else "L3F"
  comp_dir <- file.path(dm_base, comp, "Dm_color_sep")
  pops <- list(
    list(label = "Toronto", code = "HOULE"),
    list(label = "Miami",   code = "HOUSTON")
  )
  out <- map_dfr(auto_chroms, function(ch) {
    f <- file.path(comp_dir, paste0("RefAlt.", ch, ".txt"))
    if (!file.exists(f)) { cat(sprintf("  [missing] %s\n", f)); return(NULL) }
    map_dfr(pops, function(p) {
      read_line(f, paste0(p$code, "_L1F"), paste0(p$code, "_", sel_tp)) |>
        mutate(pop = p$label)
    })
  })
  apply_censor(out)
}

# ============================================================
# TEST 1: BETWEEN POPULATIONS
# ============================================================
run_between_pop <- function(dat, seed) {
  suff_arm <- make_suff(dat, c("pop", "chr"))

  # Point estimates
  pt <- map_dfr(P0_LEVELS, function(p0) {
    sig_T <- sigma_at_p0(filter(suff_arm, pop == "Toronto"), p0)
    sig_M <- sigma_at_p0(filter(suff_arm, pop == "Miami"),   p0)
    tibble(p0 = p0, sigma_toronto = sig_T, sigma_miami = sig_M,
           delta_sigma = sig_M - sig_T)
  })

  # Bootstrap: resample chromosome arms (blocks = arms)
  set.seed(seed)
  boot_d <- map_dfr(seq_len(B), function(b) {
    arms_b <- sample(auto_chroms, length(auto_chroms), replace = TRUE)
    sf_b   <- map_dfr(arms_b, function(a) filter(suff_arm, chr == a))
    map_dfr(P0_LEVELS, function(p0) {
      sig_T <- sigma_at_p0(filter(sf_b, pop == "Toronto"), p0)
      sig_M <- sigma_at_p0(filter(sf_b, pop == "Miami"),   p0)
      tibble(b = b, p0 = p0, delta = sig_M - sig_T)
    })
  })

  map_dfr(P0_LEVELS, function(p0) {
    pt_row <- filter(pt, p0 == !!p0)
    bd     <- filter(boot_d, p0 == !!p0)$delta
    bd     <- bd[is.finite(bd)]
    tibble(
      p0            = p0,
      sigma_toronto = pt_row$sigma_toronto,
      sigma_miami   = pt_row$sigma_miami,
      delta_sigma   = pt_row$delta_sigma,
      ci_lo         = quantile(bd, (1 - CONF) / 2),
      ci_hi         = quantile(bd, 1 - (1 - CONF) / 2),
      p_value       = 2 * min(mean(bd <= 0), mean(bd >= 0)),
      n_boot        = length(bd)
    )
  })
}

# ============================================================
# TEST 2: BY ARM, WITHIN POPULATION
# ============================================================
run_by_arm <- function(dat, seed) {
  # Assign 500 kb genomic blocks
  dat_blk <- dat |>
    mutate(block_id = paste0(chr, ":", (pos %/% BLOCK_SIZE) * BLOCK_SIZE))

  # Sufficient stats per (pop, chr, block_id, bin)
  suff_blk <- make_suff(dat_blk, c("pop", "chr", "block_id"))

  # Pre-split for fast bootstrap lookup: named list keyed by "pop###chr"
  suff_lookup <- split(suff_blk,
                       paste(suff_blk$pop, suff_blk$chr, sep = "###"),
                       drop = TRUE)

  get_arm_suff <- function(pop_val, chr_val)
    suff_lookup[[paste(pop_val, chr_val, sep = "###")]]

  # Resample blocks within one arm for one population
  resample_arm <- function(pop_val, chr_val) {
    sf   <- get_arm_suff(pop_val, chr_val)
    blks <- unique(sf$block_id)
    cnt  <- table(sample(blks, length(blks), replace = TRUE))
    map_dfr(names(cnt), function(bl) {
      sf[sf$block_id == bl, ] |>
        mutate(n  = n  * as.integer(cnt[[bl]]),
               s1 = s1 * as.integer(cnt[[bl]]),
               s2 = s2 * as.integer(cnt[[bl]]))
    })
  }

  # Point estimates
  pt <- map_dfr(c("Toronto", "Miami"), function(pop_val) {
    sf <- get_arm_suff(pop_val, ref_arm)
    map_dfr(P0_LEVELS, function(p0) {
      sig_ref <- sigma_at_p0(sf, p0)
      map_dfr(inverted_arms, function(arm_val) {
        sig_arm <- sigma_at_p0(get_arm_suff(pop_val, arm_val), p0)
        tibble(pop = pop_val, p0 = p0, arm = arm_val,
               sigma_arm = sig_arm, sigma_ref = sig_ref,
               delta = sig_arm - sig_ref)
      })
    })
  })

  # Bootstrap: resample blocks within each arm independently
  set.seed(seed)
  boot_d <- map_dfr(seq_len(B), function(b) {
    map_dfr(c("Toronto", "Miami"), function(pop_val) {
      # Resample each arm once, cache for all p0 levels
      resampled <- lapply(c(ref_arm, inverted_arms), function(arm_val)
        resample_arm(pop_val, arm_val))
      names(resampled) <- c(ref_arm, inverted_arms)

      map_dfr(P0_LEVELS, function(p0) {
        sig_ref <- sigma_at_p0(resampled[[ref_arm]], p0)
        map_dfr(inverted_arms, function(arm_val) {
          sig_arm <- sigma_at_p0(resampled[[arm_val]], p0)
          tibble(b = b, pop = pop_val, p0 = p0, arm = arm_val,
                 delta = sig_arm - sig_ref)
        })
      })
    })
  })

  # CIs and p-values
  map_dfr(c("Toronto", "Miami"), function(pop_val) {
    map_dfr(P0_LEVELS, function(p0) {
      map_dfr(inverted_arms, function(arm_val) {
        pt_row <- filter(pt, pop == !!pop_val, p0 == !!p0, arm == !!arm_val)
        bd     <- filter(boot_d,
                         pop == !!pop_val, p0 == !!p0, arm == !!arm_val)$delta
        bd     <- bd[is.finite(bd)]
        tibble(
          pop       = pop_val,
          p0        = p0,
          arm       = arm_val,
          sigma_arm = pt_row$sigma_arm,
          sigma_ref = pt_row$sigma_ref,
          delta     = pt_row$delta,
          ci_lo     = quantile(bd, (1 - CONF) / 2),
          ci_hi     = quantile(bd, 1 - (1 - CONF) / 2),
          p_value   = 2 * min(mean(bd <= 0), mean(bd >= 0)),
          n_boot    = length(bd)
        )
      })
    })
  })
}

# ============================================================
# MAIN
# ============================================================
bp_all  <- list()
arm_all <- list()

for (comp in comps) {
  cat(sprintf("\n=== %s ===\n", comp))
  dat <- load_dat(comp)
  n_by_pop <- dat |> count(pop)
  cat("  SNPs loaded per population:\n")
  print(as.data.frame(n_by_pop))

  cat(sprintf("  Test 1: between-population bootstrap (B = %d)...\n", B))
  bp_all[[comp]] <- run_between_pop(dat, SEED) |> mutate(comp = comp)

  cat(sprintf("  Test 2: by-arm bootstrap (B = %d)...\n", B))
  arm_all[[comp]] <- run_by_arm(dat, SEED + 1) |> mutate(comp = comp)

  rm(dat); gc()
}

bp_out  <- bind_rows(bp_all)
arm_out <- bind_rows(arm_all)

write_csv(bp_out,  file.path(out_dir, "sigma_between_pop_boot.csv"))
write_csv(arm_out, file.path(out_dir, "sigma_by_arm_boot.csv"))

# ============================================================
# PRINT TABLES
# ============================================================
fmt <- function(x, lo, hi, p) {
  sprintf("%.5f [%.5f, %.5f]  p = %.3f", x, lo, hi, p)
}

cat("\n")
cat("================================================================\n")
cat("TABLE 1: Δσ = σ(Miami) − σ(Toronto), autosomes, 90% CI\n")
cat("Prediction under inversion-inflation hypothesis: Δσ > 0\n")
cat("================================================================\n")
bp_out |>
  mutate(result = fmt(delta_sigma, ci_lo, ci_hi, p_value),
         sig    = ifelse(p_value < (1 - CONF), "*", "")) |>
  select(comp, p0, sigma_toronto, sigma_miami, delta_sigma, ci_lo, ci_hi,
         p_value, sig) |>
  arrange(comp, p0) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 5)

cat("\n")
cat("================================================================\n")
cat("TABLE 2: σ(arm) − σ(chr2R), per population, 90% CI\n")
cat("Inverted in Miami: chr2L (In(2L)t), chr3L (In(3L)P), chr3R (In(3R)K)\n")
cat("Toronto negative control: no arm carries inversions\n")
cat("================================================================\n")
arm_out |>
  mutate(sig = ifelse(p_value < (1 - CONF), "*", "")) |>
  select(comp, pop, p0, arm, sigma_arm, sigma_ref, delta, ci_lo, ci_hi,
         p_value, sig) |>
  arrange(comp, pop, p0, arm) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 5)

# ============================================================
# INTERPRETATION
# ============================================================
cat("\n================================================================\n")
cat("INTERPRETATION\n")
cat("================================================================\n")

for (comp in comps) {
  cat(sprintf("\n--- %s ---\n", gsub("_", " ", comp)))

  bp_c  <- filter(bp_out,  comp == !!comp)
  arm_c <- filter(arm_out, comp == !!comp)

  # Between-pop
  direction <- if (all(bp_c$delta_sigma < 0)) "lower" else "higher or mixed"
  sig_both  <- all(bp_c$p_value < (1 - CONF))
  max_p     <- max(bp_c$p_value)
  cat(sprintf(
    "Between populations: Miami sigma(p0) is %s than Toronto at p0 = 0.15 and 0.35",
    direction))
  if (sig_both) {
    cat(sprintf(" — significantly so at the %.0f%% level (max p = %.3f).\n",
                100 * CONF, max_p))
  } else {
    cat(sprintf(" — NOT significant at the %.0f%% level (max p = %.3f).\n",
                100 * CONF, max_p))
  }

  # Within-pop arm comparisons
  for (pop_val in c("Miami", "Toronto")) {
    arm_p <- filter(arm_c, pop == pop_val)
    n_sig <- sum(arm_p$p_value < (1 - CONF))
    cat(sprintf(
      "%s (within-arm): %d of %d arm×p0 comparisons vs chr2R significant at %.0f%% level.\n",
      pop_val, n_sig, nrow(arm_p), 100 * CONF))
  }
}

cat("\nWrote:\n")
cat("  ", file.path(out_dir, "sigma_between_pop_boot.csv"), "\n")
cat("  ", file.path(out_dir, "sigma_by_arm_boot.csv"), "\n")
