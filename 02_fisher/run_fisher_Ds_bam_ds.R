# --- GWAS on D. sim Color Exp Evo - BAM-LEVEL DOWNSAMPLED ---
# Fisher's exact test, separate tests per population (Dsim_1 and Dsim_2).
# Reads from per-comparison RefAlt files already generated from downsampled BAMs.
# No rhyper() — downsampling was done at the BAM level.
#
# Input:  process/down_sample/bam_ds/{comparison}/Ds_color_sep/RefAlt.{chr}.txt
# Output: process/down_sample/bam_ds/{comparison}/Ds_color_sep/GWAS_results_{pop}_{comp}_{chr}.txt

library(tidyverse)

# --- 1. ARGS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript run_fisher_Ds_bam_ds.R <chromosome>")
chromosome_name <- args[1]
cat("--- Chromosome:", chromosome_name, "---\n")

# --- 2. PATHS ---
refalt_base <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"

# --- 3. COMPARISON DEFINITIONS ---
comparisons <- list(
  list(name = "Dark_vs_Light",
       g1   = c("Dsim_1-L3F", "Dsim_2-L3F"),
       g2   = c("Dsim_1-L2F", "Dsim_2-L2F")),
  list(name = "Light_vs_Control",
       g1   = c("Dsim_1-L2F", "Dsim_2-L2F"),
       g2   = c("Dsim_1-L1F", "Dsim_2-L1F")),
  list(name = "Dark_vs_Control",
       g1   = c("Dsim_1-L3F", "Dsim_2-L3F"),
       g2   = c("Dsim_1-L1F", "Dsim_2-L1F"))
)

# --- 4. HELPERS ---

calc_freq <- function(ref, alt) {
  total <- ref + alt
  ifelse(total == 0, NA_real_, alt / total)
}

calc_pct_change <- function(p1, p2) {
  ifelse(is.na(p1) | is.na(p2) | p2 == 0, NA_real_, 100 * (p1 - p2) / p2)
}

run_fisher <- function(ref_g1, alt_g1, ref_g2, alt_g2) {
  n      <- length(ref_g1)
  p_vec  <- numeric(n)
  or_vec <- numeric(n)
  for (i in seq_len(n)) {
    mat <- matrix(c(ref_g1[i], alt_g1[i],
                    ref_g2[i], alt_g2[i]),
                  nrow = 2, byrow = TRUE)
    tryCatch({
      res       <- fisher.test(mat)
      p_vec[i]  <- res$p.value
      or_vec[i] <- res$estimate
    }, error = function(e) {
      p_vec[i]  <<- NA_real_
      or_vec[i] <<- NA_real_
    })
  }
  list(p = p_vec, or = or_vec)
}

# --- 5. MAIN LOOP ---
for (comp in comparisons) {
  cat("\n========================================\n")
  cat("Comparison:", comp$name, "\n")
  cat("========================================\n")

  refalt_dir  <- file.path(refalt_base, comp$name, "Ds_color_sep")
  refalt_file <- file.path(refalt_dir, paste0("RefAlt.", chromosome_name, ".txt"))

  if (!file.exists(refalt_file)) {
    cat("  Skipping — RefAlt file not found:", refalt_file, "\n")
    next
  }

  gwas_data <- read.table(refalt_file, header = TRUE, sep = "",
                           check.names = FALSE)
  cat("  Positions loaded:", nrow(gwas_data), "\n")

  # Extract REF/ALT columns for each group
  # Dsim_1 and Dsim_2 are the two populations
  sim1_g1 <- grep("Dsim_1", comp$g1, value = TRUE)
  sim1_g2 <- grep("Dsim_1", comp$g2, value = TRUE)
  sim2_g1 <- grep("Dsim_2", comp$g1, value = TRUE)
  sim2_g2 <- grep("Dsim_2", comp$g2, value = TRUE)

  sim1_ref_g1 <- gwas_data[[paste0("REF_", sim1_g1)]]
  sim1_alt_g1 <- gwas_data[[paste0("ALT_", sim1_g1)]]
  sim1_ref_g2 <- gwas_data[[paste0("REF_", sim1_g2)]]
  sim1_alt_g2 <- gwas_data[[paste0("ALT_", sim1_g2)]]

  sim2_ref_g1 <- gwas_data[[paste0("REF_", sim2_g1)]]
  sim2_alt_g1 <- gwas_data[[paste0("ALT_", sim2_g1)]]
  sim2_ref_g2 <- gwas_data[[paste0("REF_", sim2_g2)]]
  sim2_alt_g2 <- gwas_data[[paste0("ALT_", sim2_g2)]]

  cat("  Running Fisher tests...\n")
  sim1_res <- run_fisher(sim1_ref_g1, sim1_alt_g1, sim1_ref_g2, sim1_alt_g2)
  sim2_res <- run_fisher(sim2_ref_g1, sim2_alt_g1, sim2_ref_g2, sim2_alt_g2)

  # Save Dsim_1
  sim1_df <- data.frame(
    CHROM      = gwas_data$CHROM,
    POS        = gwas_data$POS,
    P_VALUE    = sim1_res$p,
    LOG10_P    = -log10(sim1_res$p),
    ODDS_RATIO = sim1_res$or,
    PCT_CHANGE = calc_pct_change(
      calc_freq(sim1_ref_g1, sim1_alt_g1),
      calc_freq(sim1_ref_g2, sim1_alt_g2)
    )
  )
  sim1_file <- file.path(refalt_dir,
    paste0("GWAS_results_Dsim_1_", comp$name, "_", chromosome_name, ".txt"))
  write.table(sim1_df, file = sim1_file, sep = " ", row.names = FALSE, quote = FALSE)
  cat("  Dsim_1 saved:", sim1_file, "\n")

  # Save Dsim_2
  sim2_df <- data.frame(
    CHROM      = gwas_data$CHROM,
    POS        = gwas_data$POS,
    P_VALUE    = sim2_res$p,
    LOG10_P    = -log10(sim2_res$p),
    ODDS_RATIO = sim2_res$or,
    PCT_CHANGE = calc_pct_change(
      calc_freq(sim2_ref_g1, sim2_alt_g1),
      calc_freq(sim2_ref_g2, sim2_alt_g2)
    )
  )
  sim2_file <- file.path(refalt_dir,
    paste0("GWAS_results_Dsim_2_", comp$name, "_", chromosome_name, ".txt"))
  write.table(sim2_df, file = sim2_file, sep = " ", row.names = FALSE, quote = FALSE)
  cat("  Dsim_2 saved:", sim2_file, "\n")

  rm(gwas_data, sim1_df, sim2_df); gc()
}

cat("\n--- Done:", chromosome_name, "---\n")
