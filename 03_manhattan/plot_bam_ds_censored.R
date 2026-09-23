# --- QQ and Manhattan plots — BAM downsampled, peak-censored ---
# Censors Dmel chr3L 8,500,000-9,000,000 and Dsim chr3L 8,100,000-8,200,000
# before any plotting. No combined plot.
# Produces per comparison:
#   - QQ plot (uncapped axes — full range shown)
#   - 5-panel Manhattan (one panel per chromosome, capped at y=100)
#   - Condensed Manhattan (all chromosomes on one x-axis, capped at y=100)

library(tidyverse)

# ============================================================
# SETTINGS
# ============================================================

comparisons <- c("Light_vs_Control", "Dark_vs_Light", "Dark_vs_Control")
chroms      <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

dm_base  <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"
ds_base  <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"
out_base <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored"

# Peak regions to remove before all plots
censor_mel <- list(chrom = "chr3L", start = 8500000, end = 9000000)
censor_sim <- list(chrom = "chr3L", start = 8342753, end = 8834090)

pop_colors <- c(
  "D. mel-1" = "dodgerblue",
  "D. mel-2" = "darkorchid2",
  "D. sim-1" = "red4",
  "D. sim-2" = "deeppink2"
)

# Italic "D. mel" / "D. sim", plain "-1" / "-2"
pop_labels <- c(
  "D. mel-1" = expression(italic("D. mel")*"-1"),
  "D. mel-2" = expression(italic("D. mel")*"-2"),
  "D. sim-1" = expression(italic("D. sim")*"-1"),
  "D. sim-2" = expression(italic("D. sim")*"-2")
)

# Candidate pigmentation genes (D. mel coordinates)
# x_nudge: horizontal offset in Mb for condensed plot only (ebony nudged left to clear burs)
candidate_genes_mel <- tibble(
  gene    = c("yellow", "tan",    "bab",    "pale",   "black",   "ebony",   "burs"),
  CHROM   = c("chrX",   "chrX",   "chr3L",  "chr3L",  "chr2L",   "chr3R",   "chr3R"),
  POS     = c(356509,   9217655,  1036369,  6713356,  13821248,  21229839,  21769518),
  x_nudge = c(0,        0,        0,        0,        0,        -1.2,       0)
)

# Candidate pigmentation genes (D. sim coordinates)
candidate_genes_sim <- tibble(
  gene    = c("yellow", "tan",    "bab",    "pale",   "black",   "ebony",   "burs"),
  CHROM   = c("chrX",   "chrX",   "chr3L",  "chr3L",  "chr2L",   "chr3R",   "chr3R"),
  POS     = c(226624,   8770795,  1044502,  6593682,  13652421,  5017916,   17995952),
  x_nudge = c(0,        0,        0,        0,        0,         0,         0)
)

# Inversions — Dmel only (dm6 coordinates)
# In(3R)K and In(3R)Payne overlap on chr3R — differentiated by fill color
inversions_mel <- tibble(
  name         = c("In(2L)t",  "In(3L)P",  "In(3R)K",  "In(3R)Payne"),
  CHROM        = c("chr2L",    "chr3L",    "chr3R",    "chr3R"),
  start        = c(2225744,    3173046,    7576895,    12257931),
  end          = c(13154180,   16302584,   25286936,   21082440),
  fill_col     = c("grey70",   "grey70",   "grey70",   "#B3CDE3"),
  y_label      = c(90,         90,         90,          82),   # condensed y
  y_label_5p   = c(96,         96,         96,          96),   # 5-panel y
  label_hjust  = c(0.5,        0.5,        0.0,         0.5)   # In(3R)K left-aligned
)

# ============================================================
# HELPERS
# ============================================================

read_gwas <- function(dir, pop, comp, chr) {
  f <- paste0(dir, "/GWAS_results_", pop, "_", comp, "_", chr, ".txt")
  if (!file.exists(f)) {
    cat(sprintf("    [missing] %s\n", f))
    return(NULL)
  }
  read_delim(f, delim = " ", show_col_types = FALSE) %>%
    mutate(CHROM = chr, POP = pop)
}

censor <- function(df, region) {
  df %>%
    filter(!(CHROM == region$chrom & POS >= region$start & POS <= region$end))
}

qq_slope <- function(log10p_vec) {
  obs <- sort(log10p_vec[!is.na(log10p_vec) & is.finite(log10p_vec) & log10p_vec >= 0],
              decreasing = TRUE)
  exp <- sort(-log10(ppoints(length(obs))), decreasing = TRUE)
  sum(obs * exp) / sum(exp^2)
}

# 5-panel Manhattan (one facet per chromosome), y capped at 100
plot_manhattan_faceted <- function(data, chrom_levels, title, ylab,
                                    gene_table = candidate_genes_mel,
                                    inv_table  = inversions_mel) {
  plot_data <- data %>%
    filter(!is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0) %>%
    mutate(CHROM   = factor(CHROM, levels = chrom_levels),
           LOG10_P = pmin(LOG10_P, 100)) %>%
    filter(!is.na(CHROM)) %>%
    slice(seq(1, n(), by = 5))

  genes_filt <- gene_table %>%
    filter(CHROM %in% chrom_levels) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_levels))

  inv_filt <- if (!is.null(inv_table) && nrow(inv_table) > 0) {
    inv_table %>%
      filter(CHROM %in% chrom_levels) %>%
      mutate(CHROM   = factor(CHROM, levels = chrom_levels),
             x_label = ifelse(label_hjust == 0, start / 1e6, (start + end) / 2 / 1e6))
  } else {
    tibble(name = character(), CHROM = factor(), start = numeric(), end = numeric(),
           fill_col = character(), y_label = numeric(), y_label_5p = numeric(),
           label_hjust = numeric(), x_label = numeric())
  }

  ggplot(plot_data, aes(x = POS / 1e6, y = LOG10_P, color = POP, alpha = LOG10_P)) +
    geom_rect(data = inv_filt,
              aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf,
                  fill = fill_col),
              alpha = 0.25, inherit.aes = FALSE) +
    scale_fill_identity() +
    geom_point(size = 0.4) +
    geom_vline(data = genes_filt,
               aes(xintercept = POS / 1e6),
               color = "grey30", linetype = "dashed", linewidth = 0.4,
               inherit.aes = FALSE) +
    geom_text(data = genes_filt,
              aes(x = POS / 1e6, y = Inf, label = gene),
              angle = 90, hjust = 1.1, vjust = 0,
              size = 3.2, color = "grey20", fontface = "italic",
              inherit.aes = FALSE) +
    geom_text(data = inv_filt,
              aes(x = x_label, y = y_label_5p, label = name, hjust = label_hjust),
              vjust = 1, size = 2.8, color = "grey30",
              fontface = "italic", inherit.aes = FALSE) +
    scale_color_manual(values = pop_colors, labels = pop_labels, name = NULL) +
    scale_alpha_continuous(range = c(0.01, 0.85), guide = "none") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    facet_grid(CHROM ~ ., scales = "free_x") +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    labs(x = "Position (Mb)", y = ylab, title = title) +
    theme_bw(base_size = 13) +
    theme(strip.text.y     = element_text(angle = 0, face = "bold"),
          legend.position  = "top",
          legend.text      = element_text(size = 13),
          axis.text        = element_text(size = 12),
          axis.title       = element_text(size = 13),
          panel.grid.minor = element_blank(),
          panel.spacing    = unit(0.5, "lines"),
          plot.margin      = margin(t = 50, r = 5, b = 5, l = 5))
}

# Condensed Manhattan (all chromosomes on one x-axis), y capped at 100
plot_manhattan_condensed <- function(data, chrom_levels, title, ylab,
                                      gene_table = candidate_genes_mel,
                                      inv_table  = inversions_mel) {

  chrom_sizes <- data %>%
    filter(!is.na(LOG10_P), is.finite(LOG10_P)) %>%
    filter(CHROM %in% chrom_levels) %>%
    group_by(CHROM) %>%
    summarise(max_pos = max(POS, na.rm = TRUE), .groups = "drop") %>%
    mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
    arrange(CHROM) %>%
    mutate(offset = lag(cumsum(as.numeric(max_pos) + 2e6), default = 0))

  plot_data <- data %>%
    filter(!is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0) %>%
    mutate(CHROM   = factor(CHROM, levels = chrom_levels),
           LOG10_P = pmin(LOG10_P, 100)) %>%
    filter(!is.na(CHROM)) %>%
    left_join(chrom_sizes %>% select(CHROM, offset), by = "CHROM") %>%
    mutate(cum_pos = POS + offset) %>%
    slice(seq(1, n(), by = 5))

  chrom_mids <- plot_data %>%
    group_by(CHROM) %>%
    summarise(mid = (min(cum_pos) + max(cum_pos)) / 2, .groups = "drop")

  genes_filt <- gene_table %>%
    filter(CHROM %in% chrom_levels) %>%
    mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
    left_join(chrom_sizes %>% select(CHROM, offset), by = "CHROM") %>%
    mutate(cum_pos = POS + offset)

  inv_filt <- if (!is.null(inv_table) && nrow(inv_table) > 0) {
    inv_table %>%
      filter(CHROM %in% chrom_levels) %>%
      mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
      left_join(chrom_sizes %>% select(CHROM, offset), by = "CHROM") %>%
      mutate(cum_start = start + offset, cum_end = end + offset,
             x_label   = ifelse(label_hjust == 0,
                                cum_start / 1e6,
                                (cum_start + cum_end) / 2 / 1e6))
  } else {
    tibble(name = character(), CHROM = factor(), cum_start = numeric(), cum_end = numeric(),
           fill_col = character(), y_label = numeric(), y_label_5p = numeric(),
           label_hjust = numeric(), x_label = numeric())
  }

  ggplot(plot_data, aes(x = cum_pos / 1e6, y = LOG10_P, color = POP, alpha = LOG10_P)) +
    geom_rect(data = inv_filt,
              aes(xmin = cum_start / 1e6, xmax = cum_end / 1e6, ymin = -Inf, ymax = Inf,
                  fill = fill_col),
              alpha = 0.25, inherit.aes = FALSE) +
    scale_fill_identity() +
    geom_point(size = 0.4) +
    geom_vline(data = genes_filt,
               aes(xintercept = cum_pos / 1e6),
               color = "grey30", linetype = "dashed", linewidth = 0.4,
               inherit.aes = FALSE) +
    geom_text(data = genes_filt,
              aes(x = cum_pos / 1e6 + x_nudge, y = Inf, label = gene),
              angle = 90, hjust = 1.1, vjust = 1,
              size = 3.2, color = "grey20", fontface = "italic",
              inherit.aes = FALSE) +
    geom_text(data = inv_filt,
              aes(x = x_label, y = y_label, label = name, hjust = label_hjust),
              vjust = 1, size = 2.8, color = "grey30",
              fontface = "italic", inherit.aes = FALSE) +
    scale_color_manual(values = pop_colors, labels = pop_labels, name = NULL) +
    scale_alpha_continuous(range = c(0.01, 0.85), guide = "none") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    scale_x_continuous(breaks = chrom_mids$mid / 1e6,
                       labels = as.character(chrom_mids$CHROM)) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    labs(x = "Chromosome", y = ylab, title = title) +
    theme_bw(base_size = 13) +
    theme(legend.position  = "top",
          legend.text      = element_text(size = 13),
          axis.text        = element_text(size = 12),
          axis.title       = element_text(size = 13),
          panel.grid.minor = element_blank(),
          plot.margin      = margin(t = 50, r = 5, b = 5, l = 5))
}

# ============================================================
# MAIN LOOP
# ============================================================

dm_pops <- c("HOULE", "HOUSTON")
ds_pops <- c("Dsim_1", "Dsim_2")

for (comp in comparisons) {
  cat("\n========================================\n")
  cat("Comparison:", comp, "\n")
  cat("========================================\n")

  dm_dir  <- file.path(dm_base, comp, "Dm_color_sep")
  ds_dir  <- file.path(ds_base, comp, "Ds_color_sep")
  out_dir <- file.path(out_base, comp)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  comp_label <- gsub("_", " ", comp)
  ylab <- expression(-log[10](p))

  # Read all data and apply censoring
  dm_raw <- map_dfr(dm_pops, function(pop) map_dfr(chroms, ~ read_gwas(dm_dir, pop, comp, .x))) %>%
    mutate(POP = recode(POP, "HOULE" = "D. mel-1", "HOUSTON" = "D. mel-2"))
  ds_raw <- map_dfr(ds_pops, function(pop) map_dfr(chroms, ~ read_gwas(ds_dir, pop, comp, .x))) %>%
    mutate(POP = recode(POP, "Dsim_1" = "D. sim-1", "Dsim_2" = "D. sim-2"))

  dm <- censor(dm_raw, censor_mel)
  ds <- censor(ds_raw, censor_sim)

  cat(sprintf("  Dmel rows after censor: %d (removed %d)\n",
              nrow(dm), nrow(dm_raw) - nrow(dm)))
  cat(sprintf("  Dsim rows after censor: %d (removed %d)\n",
              nrow(ds), nrow(ds_raw) - nrow(ds)))

  # ---- QQ plots (no axis cap — full range, Dmel and Dsim separate) ----
  cat("Building QQ plots...\n")

  make_qq_data <- function(data, pops) {
    if (length(pops) == 0) {
      cat("  [skipped] QQ — no populations in data\n")
      return(tibble())
    }
    map_dfr(pops, function(pop) {
      gwas <- data %>% filter(POP == pop, !is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0)
      if (nrow(gwas) == 0) {
        cat(sprintf("  [skipped] QQ — no valid LOG10_P rows for %s\n", pop))
        return(NULL)
      }
      s <- qq_slope(gwas$LOG10_P)
      cat(sprintf("  %s | n = %d | slope = %.4f\n", pop, nrow(gwas), s))
      obs <- sort(gwas$LOG10_P, decreasing = TRUE)
      exp <- sort(-log10(ppoints(length(obs))), decreasing = TRUE)
      data.frame(expected = exp, observed = obs, pop = pop)
    })
  }

  qq_dm <- make_qq_data(dm, unique(dm$POP))
  qq_ds <- make_qq_data(ds, unique(ds$POP))

  save_qq <- function(qq_data, label, filename) {
    if (nrow(qq_data) == 0) return(invisible(NULL))
    p <- ggplot(qq_data, aes(x = expected, y = observed, color = pop, group = pop)) +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      scale_color_manual(values = pop_colors, name = NULL) +
      guides(color = guide_legend(override.aes = list(linewidth = 2))) +
      labs(x = expression("Expected" ~ -log[10](p)),
           y = expression("Observed" ~ -log[10](p)),
           title = paste("QQ -", label, "- censored")) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    ggsave(filename, p, width = 7, height = 7, dpi = 150)
    cat(sprintf("  Saved %s\n", basename(filename)))
    rm(p); gc()
  }

  save_qq(qq_dm, paste("D. mel", comp_label),
          file.path(out_dir, paste0("QQ_Dmel_", comp, "_censored.png")))
  save_qq(qq_ds, paste("D. sim", comp_label),
          file.path(out_dir, paste0("QQ_Dsim_", comp, "_censored.png")))

  # ---- Manhattan plots ----
  cat("Building Manhattan plots...\n")

  # a. Dmel 5-panel
  if (nrow(dm) > 0) {
    p <- plot_manhattan_faceted(
      dm, chroms,
      paste0("D. mel - ", comp_label, " (censored)"), ylab,
      gene_table = candidate_genes_mel,
      inv_table  = inversions_mel
    )
    ggsave(file.path(out_dir, paste0(comp, "_Dmel_5panel_censored.png")),
           p, width = 10, height = 12, dpi = 150)
    cat("  Saved Dmel 5-panel\n")
    rm(p); gc()
  }

  # b. Dmel condensed
  if (nrow(dm) > 0) {
    p <- plot_manhattan_condensed(
      dm, chroms,
      paste0("D. mel - ", comp_label, " (censored, condensed)"), ylab,
      gene_table = candidate_genes_mel,
      inv_table  = inversions_mel
    )
    ggsave(file.path(out_dir, paste0(comp, "_Dmel_condensed_censored.png")),
           p, width = 14, height = 5, dpi = 150)
    cat("  Saved Dmel condensed\n")
    rm(p); gc()
  }

  # c. Dsim 5-panel — native Dsim coords, no inversion shading
  if (nrow(ds) > 0) {
    p <- plot_manhattan_faceted(
      ds, chroms,
      paste0("D. sim - ", comp_label, " (censored, Dsim coords)"), ylab,
      gene_table = candidate_genes_sim,
      inv_table  = NULL
    )
    ggsave(file.path(out_dir, paste0(comp, "_Dsim_5panel_censored.png")),
           p, width = 10, height = 12, dpi = 150)
    cat("  Saved Dsim 5-panel\n")
    rm(p); gc()
  }

  # d. Dsim condensed — native Dsim coords, no inversion shading
  if (nrow(ds) > 0) {
    p <- plot_manhattan_condensed(
      ds, chroms,
      paste0("D. sim - ", comp_label, " (censored, condensed, Dsim coords)"), ylab,
      gene_table = candidate_genes_sim,
      inv_table  = NULL
    )
    ggsave(file.path(out_dir, paste0(comp, "_Dsim_condensed_censored.png")),
           p, width = 14, height = 5, dpi = 150)
    cat("  Saved Dsim condensed\n")
    rm(p); gc()
  }

  rm(dm, ds, dm_raw, ds_raw); gc()
  cat("Done:", comp, "\n")
}

cat("\n--- All done ---\n")
