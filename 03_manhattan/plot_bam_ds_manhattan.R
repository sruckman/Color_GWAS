# --- QQ and Manhattan plots for BAM-LEVEL DOWNSAMPLED Fisher results ---
# Reads from process/down_sample/bam_ds/{comparison}/Dm_color_sep/ and Ds_color_sep/
# No p-value correction — downsampling replaces that step.
# Produces per-comparison QQ plots and Manhattan plots
# (Dmel only, Dsim only, combined with gene annotations and inversion shading).

library(tidyverse)

# ============================================================
# SETTINGS
# ============================================================

comparisons <- c("Light_vs_Control", "Dark_vs_Light", "Dark_vs_Control")
chroms      <- c("chr2L", "chr2R", "chr3L", "chr3R", "chrX")

dm_base <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"
ds_base <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"
out_base <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots"

pop_colors <- c(
  "HOULE"   = "dodgerblue",
  "HOUSTON" = "darkorchid2",
  "Dsim_1"  = "red4",
  "Dsim_2"  = "deeppink2"
)

pop_labels <- c(
  HOULE   = "Toronto",
  HOUSTON = "Miami",
  Dsim_1  = "Tallahassee-1",
  Dsim_2  = "Tallahassee-2"
)

# Candidate pigmentation genes (D. mel coordinates)
# hjust_label: controls label height in the top margin for angle=90 text.
# Stagger ebony and burs (540 kb apart on chr3R) so labels don't overlap.
candidate_genes_mel <- tibble(
  gene        = c("yellow", "tan",    "bab",    "pale",   "black",   "ebony",   "burs"),
  CHROM       = c("chrX",   "chrX",   "chr3L",  "chr3L",  "chr2L",   "chr3R",   "chr3R"),
  POS         = c(356509,   9217655,  1036369,  6713356,  13821248,  21229839,  21769518),
  hjust_label = c(1.1,      1.1,      1.1,      1.1,      1.1,       1.1,       1.7)
)

# Candidate pigmentation genes (D. sim coordinates)
candidate_genes_sim <- tibble(
  gene        = c("yellow", "tan",    "bab",    "pale",   "black",   "ebony",   "burs"),
  CHROM       = c("chrX",   "chrX",   "chr3L",  "chr3L",  "chr2L",   "chr3R",   "chr3R"),
  POS         = c(226624,   8770795,  1044502,  6593682,  13652421,  5017916,   17995952),
  hjust_label = c(1.1,      1.1,      1.1,      1.1,      1.1,       1.1,       1.1)
)

# Polymorphic inversions in D. mel only
# In(2L)t: chr2L:2,225,744-13,154,180  (dm6, confirmed FlyBase FBab0004696)
# In(3R)P: chr3R:12,257,931-21,082,440 (dm6, Corbett-Detig & Hartl 2012)
inversions_mel <- tibble(
  name       = c("In(2L)t",  "In(3L)P",  "In(3R)K",  "In(3R)P"),
  CHROM      = c("chr2L",    "chr3L",     "chr3R",     "chr3R"),
  start      = c(2225744L,   3173046L,    7576289L,    12257931L),
  end        = c(13154180L,  16301941L,   22422742L,   21082440L),
  fill_col   = c("grey70",   "grey70",    "grey70",    "#ADD8E6"),  # In(3R)P = light blue
  label_side = c("center",   "right",     "left",      "center")   # In(3L)P right, In(3R)K left
)

# Chromosome sizes for cumulative-position Manhattan layout.
# Dmel: dm6 lengths. Dsim: max sim_end values from sim_to_mel_map.
chrom_sizes_mel <- c(
  chrX   = 23542271L, chr2L = 23513712L, chr2R = 25286936L,
  chr3L  = 28110227L, chr3R = 32079331L
)
chrom_sizes_sim <- c(
  chrX   = 22032349L, chr2L = 23642019L, chr2R = 22319025L,
  chr3L  = 23379930L, chr3R = 28132288L
)

# Chorion gene cluster windows — NATIVE COORDINATES per species.
# Flagged during data read so the flag survives map_dsim().
# D. mel:  3L 66D (8.5-9.0 Mb) and X 7F (8.4-8.6 Mb; Cp36 at chrX:8,479,319-8,480,523 dm6)
# D. sim:  3L 66D (8.34-8.83 Mb) and X 7F (8.70-8.78 Mb; from coverage scan)
# IMPORTANT: D. sim tan is at chrX 8,770,795. The X window ends at 8,784,000 to stay
# below tan + a small buffer. Do NOT widen the Dsim X window without re-checking tan.
chorion_mel <- tribble(
  ~CHROM,  ~start,    ~end,
  "chr3L", 8500000L,  9000000L,
  "chrX",  8400000L,  8600000L
)
chorion_sim <- tribble(
  ~CHROM,  ~start,    ~end,
  "chr3L", 8342753L,  8834090L,
  "chrX",  7800000L,  8400000L   # Confirmed by coverage scan (chrX bracket 7.0-8.5 Mb).
                                  # Peak at 8,208,001 in 5/6 lines (up to 68x amplified);
                                  # secondary peak at 7,860,001 in Dsim_1-L3F.
                                  # Old coordinates (8,702-8,784 kb) were near tan (8,770,795)
                                  # and missed the cluster entirely.
)

# Dsim -> Dmel coordinate mapping
sim_to_mel_map <- tibble(
  sim_chr   = c("chr2L",  "chr2R",  "chr3L", "chr3L",  "chr3L",  "chr3L",  "chr3R",  "chr3R",  "chr3R",  "chr3R",  "chr3R",  "chrX"),
  sim_start = c( 78303,    29807,    31826,   7851416,  8834091,  19666517,  129992,  5252305,  18227079,  4881699,  26119855,  33691),
  sim_end   = c(23642019, 22319025, 7851415,  8834090, 19666516,  23379930, 3487887, 17970306,  26119854,  5252304,  28132288, 22032349),
  mel_chr   = c("chr2L",  "chr2R",  "chr3L", "chr3L",  "chr3L",  "chr3L",  "chr3R",  "chr3R",  "chr3R",  "chr3R",  "chr3R",  "chrX"),
  mel_start = c(       0,        0,       0, 8000000,  9000001,  20000001,        0, 7000001,  22000000,  21000001, 30000001,       0),
  mel_end   = c(25000000, 26000000, 7999999, 9000000, 20000000,  25000000,  7000000, 21000000, 30000000,  22000000, 32500000, 25000000),
  inverted  = c(FALSE,    FALSE,    FALSE,   FALSE,    FALSE,    FALSE,    FALSE,    FALSE,    FALSE,    TRUE,     FALSE,    FALSE)
)

# ============================================================
# HELPERS
# ============================================================

read_gwas <- function(dir, pop, comp, chr) {
  f <- paste0(dir, "/GWAS_results_", pop, "_", comp, "_", chr, ".txt")
  if (!file.exists(f)) return(NULL)
  read_delim(f, delim = " ", show_col_types = FALSE) %>%
    mutate(CHROM = chr, POP = pop)
}

# Flag SNPs overlapping chorion amplification windows in native coordinates.
# Must be called before map_dsim() so Dsim points use Dsim coordinates.
# chorion_tbl has columns: CHROM, start, end.
flag_chorion <- function(df, chorion_tbl) {
  is_ch <- rep(FALSE, nrow(df))
  for (i in seq_len(nrow(chorion_tbl))) {
    is_ch <- is_ch | (df$CHROM == chorion_tbl$CHROM[i] &
                        df$POS  >= chorion_tbl$start[i] &
                        df$POS  <= chorion_tbl$end[i])
  }
  df %>% mutate(is_chorion = is_ch)
}

# Map Dsim native coordinates to Dmel coordinates.
# Carries is_chorion so greying survives the remapping.
map_dsim <- function(df) {
  df %>%
    rename(sim_chr = CHROM) %>%
    left_join(sim_to_mel_map, by = "sim_chr") %>%
    filter(POS >= sim_start, POS <= sim_end) %>%
    mutate(
      CHROM = mel_chr,
      POS   = if_else(
        inverted,
        mel_end   - (POS - sim_start) * ((mel_end - mel_start) / (sim_end - sim_start)),
        mel_start + (POS - sim_start) * ((mel_end - mel_start) / (sim_end - sim_start))
      )
    ) %>%
    select(CHROM, POS, LOG10_P, POP, is_chorion)
}

plot_manhattan <- function(data, chrom_levels, title, ylab,
                            gene_table  = candidate_genes_mel,
                            inv_table   = inversions_mel,
                            chrom_sizes = chrom_sizes_mel) {
  GAP <- 2e6L

  plot_data <- data %>%
    filter(!is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0) %>%
    filter(CHROM %in% chrom_levels) %>%
    group_by(CHROM, POP) %>%
    arrange(POS, .by_group = TRUE) %>%
    slice(seq(1L, n(), by = 5L)) %>%
    ungroup()

  # is_chorion flagged in native coords before this function — do not recompute.
  if (!"is_chorion" %in% names(plot_data))
    plot_data <- plot_data %>% mutate(is_chorion = FALSE)

  # Cumulative positions
  sizes   <- chrom_sizes[chrom_levels]
  offsets <- cumsum(c(0L, head(sizes + GAP, -1L)))
  names(offsets) <- chrom_levels
  plot_data <- plot_data %>% mutate(cum_pos = POS + offsets[CHROM])

  # Thin vertical lines at chromosome boundaries (between arms)
  chrom_bounds <- tibble(x = (offsets[-1] - GAP / 2L) / 1e6)

  # X-axis: chromosome name at midpoint
  chrom_mids   <- (offsets + sizes[chrom_levels] / 2) / 1e6
  chrom_labels <- gsub("chr", "", chrom_levels)

  genes_filt <- gene_table %>%
    filter(CHROM %in% chrom_levels) %>%
    mutate(cum_pos = POS + offsets[CHROM])

  inv_filt <- if (!is.null(inv_table)) {
    inv_table %>%
      filter(CHROM %in% chrom_levels) %>%
      mutate(
        cum_start = start + offsets[CHROM],
        cum_end   = end   + offsets[CHROM],
        # Per-inversion label x position and justification
        label_x   = case_when(
          label_side == "left"  ~ cum_start / 1e6,
          label_side == "right" ~ cum_end   / 1e6,
          TRUE                  ~ (cum_start + cum_end) / 2 / 1e6
        ),
        label_hjust = case_when(
          label_side == "left"  ~ 0,
          label_side == "right" ~ 1,
          TRUE                  ~ 0.5
        )
      )
  } else {
    tibble(name = character(), CHROM = character(),
           cum_start = numeric(), cum_end = numeric(),
           fill_col = character(), label_x = numeric(), label_hjust = numeric())
  }

  inv_grey <- filter(inv_filt, fill_col == "grey70")
  inv_blue <- filter(inv_filt, fill_col != "grey70")

  ggplot(plot_data, aes(x = cum_pos / 1e6, y = LOG10_P, color = POP)) +
    # Chromosome boundary lines (subtle)
    geom_vline(data = chrom_bounds, aes(xintercept = x),
               colour = "grey80", linewidth = 0.4, inherit.aes = FALSE) +
    # Inversion shading — grey and light-blue inversions drawn separately
    geom_rect(data = inv_grey,
              aes(xmin = cum_start / 1e6, xmax = cum_end / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "grey70", alpha = 0.18, inherit.aes = FALSE) +
    geom_rect(data = inv_blue,
              aes(xmin = cum_start / 1e6, xmax = cum_end / 1e6,
                  ymin = -Inf, ymax = Inf),
              fill = "#ADD8E6", alpha = 0.35, inherit.aes = FALSE) +
    # Inversion labels — white background, per-inversion x position and hjust
    geom_label(data = inv_filt,
               aes(x = label_x, y = Inf, label = name, hjust = label_hjust),
               vjust = 1.3, size = 3.2,
               colour = "grey25", fill = "white", label.size = 0,
               label.padding = unit(0.12, "lines"), inherit.aes = FALSE) +
    # All colored points first (non-chorion)
    geom_point(data = dplyr::filter(plot_data, !is_chorion),
               size = 0.5, alpha = 0.45) +
    # Chorion grey points plotted LAST so they sit on top and pop out
    geom_point(data = dplyr::filter(plot_data, is_chorion),
               aes(x = cum_pos / 1e6, y = LOG10_P),
               size = 0.5, colour = "grey45", alpha = 0.85, inherit.aes = FALSE) +
    # Gene dashed lines
    geom_vline(data = genes_filt,
               aes(xintercept = cum_pos / 1e6),
               colour = "grey30", linetype = "dashed", linewidth = 0.4,
               inherit.aes = FALSE) +
    # Gene labels — white fill box covers the line.
    # Three layers: standard genes centered, ebony left of line, burs right of line.
    # vjust with angle=90: 0.5=centered, >0.5 shifts LEFT, <0.5 shifts RIGHT.
    geom_label(data = filter(genes_filt, !gene %in% c("ebony", "burs")),
               aes(x = cum_pos / 1e6, y = Inf, label = gene),
               angle = 90, hjust = 1.1, vjust = 0.5,
               size = 3.0, colour = "grey20", fill = "white", label.size = 0,
               label.padding = unit(0.08, "lines"), inherit.aes = FALSE) +
    # ebony: label to the LEFT of its dashed line
    geom_label(data = filter(genes_filt, gene == "ebony"),
               aes(x = cum_pos / 1e6, y = Inf, label = gene),
               angle = 90, hjust = 1.1, vjust = 1.1,
               size = 3.0, colour = "grey20", fill = "white", label.size = 0,
               label.padding = unit(0.08, "lines"), inherit.aes = FALSE) +
    # burs: label to the RIGHT of its dashed line (~540 kb from ebony on chr3R)
    geom_label(data = filter(genes_filt, gene == "burs"),
               aes(x = cum_pos / 1e6, y = Inf, label = gene),
               angle = 90, hjust = 1.1, vjust = -0.1,
               size = 3.0, colour = "grey20", fill = "white", label.size = 0,
               label.padding = unit(0.08, "lines"), inherit.aes = FALSE) +
    scale_x_continuous(breaks = chrom_mids, labels = chrom_labels,
                        expand = expansion(mult = 0.01)) +
    scale_color_manual(values = pop_colors, name = NULL, labels = pop_labels) +
    guides(colour = guide_legend(
      override.aes = list(size = 5, alpha = 1)   # large visible dots in legend
    )) +
    coord_cartesian(clip = "off") +
    labs(x = "Chromosome", y = ylab, title = title) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "top",
      legend.text      = element_text(size = 13),
      legend.key.size  = unit(1.0, "cm"),
      axis.text        = element_text(size = 12),
      axis.title       = element_text(size = 13),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "grey75"),
      plot.margin        = margin(t = 30, r = 5, b = 5, l = 5)
    )
}

qq_slope <- function(log10p_vec) {
  obs <- sort(log10p_vec[!is.na(log10p_vec) & is.finite(log10p_vec) & log10p_vec >= 0],
              decreasing = TRUE)
  exp <- sort(-log10(ppoints(length(obs))), decreasing = TRUE)
  sum(obs * exp) / sum(exp^2)
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

  # ---- QQ plot ----
  cat("Building QQ plot...\n")
  qq_data <- map_dfr(c(
    lapply(dm_pops, function(pop) {
      gwas <- map_dfr(chroms, ~ read_gwas(dm_dir, pop, comp, .x)) %>%
        filter(!is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0)
      if (nrow(gwas) == 0) return(NULL)
      s <- qq_slope(gwas$LOG10_P)
      cat(sprintf("  %s | slope = %.4f\n", pop, s))
      obs <- sort(gwas$LOG10_P, decreasing = TRUE)
      exp <- sort(-log10(ppoints(length(obs))), decreasing = TRUE)
      data.frame(expected = exp, observed = obs, pop = pop) %>%
        mutate(row = row_number()) %>%
        filter(observed >= 5 | row %% 50 == 0)
    }),
    lapply(ds_pops, function(pop) {
      gwas <- map_dfr(chroms, ~ read_gwas(ds_dir, pop, comp, .x)) %>%
        filter(!is.na(LOG10_P), is.finite(LOG10_P), LOG10_P >= 0)
      if (nrow(gwas) == 0) return(NULL)
      s <- qq_slope(gwas$LOG10_P)
      cat(sprintf("  %s | slope = %.4f\n", pop, s))
      obs <- sort(gwas$LOG10_P, decreasing = TRUE)
      exp <- sort(-log10(ppoints(length(obs))), decreasing = TRUE)
      data.frame(expected = exp, observed = obs, pop = pop) %>%
        mutate(row = row_number()) %>%
        filter(observed >= 5 | row %% 50 == 0)
    })
  ), identity)

  if (nrow(qq_data) > 0) {
    p <- ggplot(qq_data, aes(x = expected, y = observed, color = pop, group = pop)) +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      scale_color_manual(values = pop_colors, name = NULL, labels = pop_labels) +
      labs(x = expression("Expected" ~ -log[10](p)),
           y = expression("Observed" ~ -log[10](p)),
           title = paste("QQ -", comp_label, "- BAM downsampled")) +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    ggsave(file.path(out_dir, paste0("QQ_", comp, ".png")),
           p, width = 7, height = 7, dpi = 150)
    cat("  Saved QQ\n")
    rm(p); gc()
  }

  # ---- Manhattan plots ----
  cat("Building Manhattan plots...\n")

  dm <- map_dfr(dm_pops, function(pop) map_dfr(chroms, ~ read_gwas(dm_dir, pop, comp, .x)))
  ds <- map_dfr(ds_pops, function(pop) map_dfr(chroms, ~ read_gwas(ds_dir, pop, comp, .x)))

  # Flag chorion windows in each species' native coordinates, before any coord remapping.
  if (nrow(dm) > 0) dm <- flag_chorion(dm, chorion_mel)
  if (nrow(ds) > 0) ds <- flag_chorion(ds, chorion_sim)

  cat(sprintf("  Dmel rows: %d  Dsim rows: %d\n", nrow(dm), nrow(ds)))

  # a. Dmel only
  if (nrow(dm) > 0) {
    p <- plot_manhattan(dm, chroms,
           paste0("D. mel - ", comp_label, " (BAM downsampled)"), ylab)
    ggsave(file.path(out_dir, paste0(comp, "_Dmel.png")),
           p, width = 16, height = 5, dpi = 150)
    cat("  Saved Dmel Manhattan\n")
    rm(p); gc()
  }

  # b. Dsim only — native Dsim coords, Dsim gene table, no inversion shading
  if (nrow(ds) > 0) {
    p <- plot_manhattan(ds, chroms,
           paste0("D. sim - ", comp_label, " (BAM downsampled, Dsim coords)"), ylab,
           gene_table  = candidate_genes_sim,
           inv_table   = NULL,
           chrom_sizes = chrom_sizes_sim)
    ggsave(file.path(out_dir, paste0(comp, "_Dsim.png")),
           p, width = 16, height = 5, dpi = 150)
    cat("  Saved Dsim Manhattan\n")
    rm(p); gc()
  }

  # c. Combined — Dsim mapped to mel coords
  if (nrow(dm) > 0 && nrow(ds) > 0) {
    ds_mapped <- map_dsim(ds)
    combined  <- bind_rows(dm, ds_mapped)
    p <- plot_manhattan(combined, chroms,
           paste0("Combined - ", comp_label, " (BAM downsampled)"), ylab)
    ggsave(file.path(out_dir, paste0(comp, "_combined.png")),
           p, width = 16, height = 5, dpi = 150)
    cat("  Saved combined Manhattan\n")
    rm(p, ds_mapped, combined); gc()
  }

  rm(dm, ds); gc()
}

cat("\nDone.\n")
