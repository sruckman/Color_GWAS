# --- GWAS on D. mel Color Exp Evo - BAM-LEVEL DOWNSAMPLED ---
# Fisher's exact test, separate tests per population (HOULE and HOUSTON).
# Reads from per-comparison RefAlt files already generated from downsampled BAMs.
# No rhyper() — downsampling was done at the BAM level.
#
# Input:  process/down_sample/bam_ds/{comparison}/Dm_color_sep/RefAlt.{chr}.txt
# Output: process/down_sample/bam_ds/{comparison}/Dm_color_sep/GWAS_results_{pop}_{comp}_{chr}.txt

library(tidyverse)

# --- 1. ARGS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Usage: Rscript run_fisher_Dm_bam_ds.R <chromosome>")
chromosome_name <- args[1]
cat("--- Chromosome:", chromosome_name, "---\n")

# --- 2. PATHS ---
refalt_base <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"

# --- 3. COMPARISON DEFINITIONS ---
comparisons <- list(
  list(name = "Dark_vs_Light",
       g1   = c("HOULE_L3F",  "HOUSTON_L3F"),
       g2   = c("HOULE_L2F",  "HOUSTON_L2F")),
  list(name = "Light_vs_Control",
       g1   = c("HOULE_L2F",  "HOUSTON_L2F"),
       g2   = c("HOULE_L1F",  "HOUSTON_L1F")),
  list(name = "Dark_vs_Control",
       g1   = c("HOULE_L3F",  "HOUSTON_L3F"),
       g2   = c("HOULE_L1F",  "HOUSTON_L1F"))
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

  refalt_dir  <- file.path(refalt_base, comp$name, "Dm_color_sep")
  refalt_file <- file.path(refalt_dir, paste0("RefAlt.", chromosome_name, ".txt"))

  if (!file.exists(refalt_file)) {
    cat("  Skipping — RefAlt file not found:", refalt_file, "\n")
    next
  }

  gwas_data <- read.table(refalt_file, header = TRUE, sep = "",
                           check.names = FALSE)
  cat("  Positions loaded:", nrow(gwas_data), "\n")

  # Extract REF/ALT columns for each group
  houle_g1   <- grep("HOULE",   comp$g1, value = TRUE)
  houle_g2   <- grep("HOULE",   comp$g2, value = TRUE)
  houston_g1 <- grep("HOUSTON", comp$g1, value = TRUE)
  houston_g2 <- grep("HOUSTON", comp$g2, value = TRUE)

  houle_ref_g1   <- gwas_data[[paste0("REF_", houle_g1)]]
  houle_alt_g1   <- gwas_data[[paste0("ALT_", houle_g1)]]
  houle_ref_g2   <- gwas_data[[paste0("REF_", houle_g2)]]
  houle_alt_g2   <- gwas_data[[paste0("ALT_", houle_g2)]]

  houston_ref_g1 <- gwas_data[[paste0("REF_", houston_g1)]]
  houston_alt_g1 <- gwas_data[[paste0("ALT_", houston_g1)]]
  houston_ref_g2 <- gwas_data[[paste0("REF_", houston_g2)]]
  houston_alt_g2 <- gwas_data[[paste0("ALT_", houston_g2)]]

  cat("  Running Fisher tests...\n")
  houle_res   <- run_fisher(houle_ref_g1,   houle_alt_g1,   houle_ref_g2,   houle_alt_g2)
  houston_res <- run_fisher(houston_ref_g1, houston_alt_g1, houston_ref_g2, houston_alt_g2)

  # Save HOULE
  houle_df <- data.frame(
    CHROM      = gwas_data$CHROM,
    POS        = gwas_data$POS,
    P_VALUE    = houle_res$p,
    LOG10_P    = -log10(houle_res$p),
    ODDS_RATIO = houle_res$or,
    PCT_CHANGE = calc_pct_change(
      calc_freq(houle_ref_g1, houle_alt_g1),
      calc_freq(houle_ref_g2, houle_alt_g2)
    )
  )
  houle_file <- file.path(refalt_dir,
    paste0("GWAS_results_HOULE_", comp$name, "_", chromosome_name, ".txt"))
  write.table(houle_df, file = houle_file, sep = " ", row.names = FALSE, quote = FALSE)
  cat("  HOULE saved:", houle_file, "\n")

  # Save HOUSTON
  houston_df <- data.frame(
    CHROM      = gwas_data$CHROM,
    POS        = gwas_data$POS,
    P_VALUE    = houston_res$p,
    LOG10_P    = -log10(houston_res$p),
    ODDS_RATIO = houston_res$or,
    PCT_CHANGE = calc_pct_change(
      calc_freq(houston_ref_g1, houston_alt_g1),
      calc_freq(houston_ref_g2, houston_alt_g2)
    )
  )
  houston_file <- file.path(refalt_dir,
    paste0("GWAS_results_HOUSTON_", comp$name, "_", chromosome_name, ".txt"))
  write.table(houston_df, file = houston_file, sep = " ", row.names = FALSE, quote = FALSE)
  cat("  HOUSTON saved:", houston_file, "\n")

  rm(gwas_data, houle_df, houston_df); gc()
}

cat("\n--- Done:", chromosome_name, "---\n")
