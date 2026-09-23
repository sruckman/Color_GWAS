#!/usr/bin/env Rscript
# run_inversion_noise_Dm.R
#
# Tony's "accidental inversion experiment".
#
# Miami (internal code HOUSTON) segregates for In(2L)t, In(3L)P, In(3R)K.
# Toronto (internal code HOULE) carries none of them.
# Both were selected the same way, in the same species, at the same time.
# So we already have an inversion vs inversion-free contrast, by accident.
#
# Question: do inversions inflate the genome-wide noise in allele frequency
# change, and therefore the resolution ceiling?
#
# Test 1 (between populations): sigma(p0) in Miami vs Toronto.
# Test 2 (within Miami):     sigma(p0) on inverted arms (2L, 3L, 3R)
#                               vs non-inverted arms (2R, X).
#                               Toronto is the negative control: no arm is
#                               inverted, so no arm contrast should appear.
#
# If inversions drive the smear, Miami should be noisier than Toronto,
# the effect should be concentrated on the inverted arms, and Toronto should
# show no arm effect.
#
# Output: process/down_sample/bam_ds/plots_censored/inversion_noise/
#   fig_sigma_by_pop.png       sigma(p0), Toronto vs Miami
#   fig_sigma_by_arm.png       sigma(p0) by arm, faceted by population
#   sigma_by_stratum.csv       sigma at p0 = 0.15 and 0.35 for every chr x pop
#   power_by_pop.csv           single-population (R=1) power, each population

suppressMessages(library(tidyverse))

# ============================================================
# PARAMETERS  (match run_design_power_Dm.R)
# ============================================================
params <- list(
  sel_frac  = c(0.800, 0.732, 0.644, 0.596, 0.528, 0.460, 0.392, 0.324,
                0.256, 0.188, 0.168, 0.148, 0.120, 0.100, 0.084, 0.084),
  n_scored  = 500,
  p0_levels = c(0.05, 0.15, 0.35),
  e_grid    = seq(0.05, 1.5, 0.025),
  n_snps    = 2e6,
  alpha     = 0.05,
  freq_binw = 0.02
)

dm_base <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"
out_dir <- file.path(dm_base, "plots_censored", "inversion_noise")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chroms <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")
censor <- list(chrom = "chr3L", start = 8500000L, end = 9000000L)

# Arms carrying a cosmopolitan inversion that is polymorphic in D. mel-2.
# In(2L)t on 2L, In(3L)P on 3L, In(3R)K and In(3R)Payne on 3R.
inverted_arms <- c("chr2L", "chr3L", "chr3R")

comps <- c("Dark_vs_Control", "Light_vs_Control")

# ============================================================
# CORE FUNCTIONS  (identical to run_design_power_Dm.R)
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
    out <- out |> filter(!(chr == censor$chrom &
                             pos >= censor$start & pos <= censor$end))
  out
}

er_signal <- function(e, p0, er_i) {
  p <- p0
  for (ig in er_i) p <- min(max(p + ig * e * p * (1 - p), 0), 1)
  p - p0
}

make_power_fn <- function(n_snps, alpha) {
  zc <- qnorm(1 - 0.5 * (alpha / n_snps))
  function(ncp) pnorm(ncp - zc) + pnorm(-ncp - zc)
}

# sigma(p0) for an arbitrary subset of SNPs
sigma_curve <- function(dat, bw = 0.02, min_n = 100) {
  mids <- seq(bw / 2, 1 - bw / 2, bw)
  dat |>
    mutate(bin = pmin(floor(p0 / bw) + 1L, length(mids))) |>
    group_by(bin) |>
    summarise(n = n(), s1 = sum(dp), s2 = sum(dp^2), .groups = "drop") |>
    filter(n >= min_n) |>
    mutate(p_mid = mids[bin],
           sd    = sqrt(pmax((s2 - s1^2 / n) / (n - 1), 0))) |>
    select(p_mid, n, sd)
}

sigma_at <- function(curve, p0) approx(curve$p_mid, curve$sd, xout = p0, rule = 2)$y

load_all_strata_Dm <- function(comp_name, chroms, censor = NULL) {
  sel_tp   <- if (grepl("Light", comp_name)) "L2F" else "L3F"
  comp_dir <- file.path(dm_base, comp_name, "Dm_color_sep")

  pops <- list(
    list(label = "Toronto", code = "HOULE"),    # high latitude, inversion-free
    list(label = "Miami",   code = "HOUSTON")   # low latitude, carries In(2L)t, In(3L)P, In(3R)K
  )

  map_dfr(chroms, function(ch) {
    f <- file.path(comp_dir, paste0("RefAlt.", ch, ".txt"))
    if (!file.exists(f)) { cat(sprintf("  [missing] %s\n", basename(f))); return(NULL) }
    map_dfr(pops, function(p) {
      read_line(f, paste0(p$code, "_L1F"), paste0(p$code, "_", sel_tp), censor) |>
        mutate(pop = p$label, chr = ch)
    })
  }) |>
    mutate(arm_type = if_else(chr %in% inverted_arms,
                              "Inverted arm (2L, 3L, 3R)",
                              "Non-inverted arm (2R, X)"))
}

# ============================================================
# RUN
# ============================================================
er_i <- sel_intensity(round(params$sel_frac * params$n_scored) / params$n_scored)
powr <- make_power_fn(params$n_snps, params$alpha)

sigma_tbl <- list()
power_tbl <- list()
bypop_all <- list()
byarm_all <- list()

for (comp in comps) {
  cat("=== ", comp, " ===\n")
  dat <- load_all_strata_Dm(comp, chroms, censor)
  if (nrow(dat) == 0) next

  # ---- sigma(p0) by population ----
  bypop <- dat |>
    group_nest(pop) |>
    mutate(curve = map(data, sigma_curve)) |>
    select(-data) |>
    unnest(curve) |>
    mutate(comp = comp)
  bypop_all[[comp]] <- bypop

  # ---- sigma(p0) by population x arm type ----
  byarm <- dat |>
    group_nest(pop, arm_type) |>
    mutate(curve = map(data, sigma_curve)) |>
    select(-data) |>
    unnest(curve) |>
    mutate(comp = comp)
  byarm_all[[comp]] <- byarm

  # ---- sigma at each chromosome x population (the table Tony will want) ----
  sigma_tbl[[comp]] <- dat |>
    group_nest(pop, chr, arm_type) |>
    mutate(curve = map(data, sigma_curve),
           n_snp = map_int(data, nrow),
           sd_p15 = map_dbl(curve, sigma_at, 0.15),
           sd_p35 = map_dbl(curve, sigma_at, 0.35)) |>
    select(pop, chr, arm_type, n_snp, sd_p15, sd_p35) |>
    mutate(comp = comp)

  # ---- single-population power (R = 1), each population on its own ----
  power_tbl[[comp]] <- bypop |>
    group_nest(pop) |>
    mutate(pw = map(data, function(cv) {
      expand_grid(p0 = params$p0_levels, e = params$e_grid) |>
        mutate(sig   = sigma_at(cv, p0),
               dp_er = map2_dbl(e, p0, \(e, p) er_signal(e, p, er_i)),
               power = powr(dp_er / sig))
    })) |>
    select(-data) |>
    unnest(pw) |>
    mutate(comp = comp)
}

sigma_out <- bind_rows(sigma_tbl)
power_out <- bind_rows(power_tbl)
write_csv(sigma_out, file.path(out_dir, "sigma_by_stratum.csv"))
write_csv(power_out, file.path(out_dir, "power_by_pop.csv"))

# ---- printed summary: the actual answer to Tony's question ----
cat("\n########## sigma(p0) by population ##########\n")
sigma_out |>
  group_by(comp, pop) |>
  summarise(sd_p15 = mean(sd_p15), sd_p35 = mean(sd_p35), .groups = "drop") |>
  as.data.frame() |> print()

cat("\n########## sigma(p0) by population x arm type ##########\n")
sigma_out |>
  group_by(comp, pop, arm_type) |>
  summarise(sd_p15 = mean(sd_p15), sd_p35 = mean(sd_p35), .groups = "drop") |>
  as.data.frame() |> print()

cat("\n########## R=1 power ceiling, each population alone ##########\n")
power_out |>
  group_by(comp, pop, p0) |>
  summarise(max_power = max(power), .groups = "drop") |>
  as.data.frame() |> print()

# ============================================================
# PLOTS
# ============================================================
pop_cols <- c("Toronto" = "dodgerblue3", "Miami" = "darkorchid3")

p1 <- bind_rows(bypop_all) |>
  ggplot(aes(p_mid, sd, colour = pop)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~comp) +
  scale_colour_manual(values = pop_cols, name = NULL) +
  labs(x = "Starting (control) allele frequency",
       y = "SD of allele-frequency change",
       title = "Inversion-carrying vs inversion-free D. melanogaster populations",
       subtitle = "Miami (low latitude) segregates for In(2L)t, In(3L)P and In(3R)K.  Toronto (high latitude) carries none.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_sigma_by_pop.png"), p1,
       width = 9, height = 4.5, dpi = 300)

p2 <- bind_rows(byarm_all) |>
  ggplot(aes(p_mid, sd, colour = arm_type)) +
  geom_line(linewidth = 0.9) +
  facet_grid(comp ~ pop) +
  scale_colour_manual(values = c("Inverted arm (2L, 3L, 3R)" = "firebrick",
                                 "Non-inverted arm (2R, X)"  = "grey40"),
                      name = NULL) +
  labs(x = "Starting (control) allele frequency",
       y = "SD of allele-frequency change",
       title = "Does the noise sit on the inverted arms?",
       subtitle = "Toronto is the negative control: it carries no inversions, so no arm contrast is expected.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_sigma_by_arm.png"), p2,
       width = 9, height = 6.5, dpi = 300)

cat("\nWrote output to: ", out_dir, "\n")
