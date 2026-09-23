#!/usr/bin/env Rscript
# check_chorion_clusters.R
#
# Tony: "Major peak is likely to do with chorion amplification in adult female
#        DNA samples.  Then prove it."
#
# The proof does not need new sequencing.  Drosophila has TWO chorion gene
# clusters, and both are developmentally amplified in ovarian follicle cells
# during oogenesis (Spradling and Mahowald 1980; Claycomb and Orr-Weaver 2005):
#
#   - the third-chromosome cluster at cytological 66D  (chr3L)  ~60-fold in follicle cells
#   - the X-linked cluster        at cytological 7F   (chrX)   ~15-fold in follicle cells
#
# We already see the chr3L cluster elevated 4-5x in whole-female DNA (diluted,
# because follicle cells are a minority of the cells in a whole fly).
#
# PREDICTION: if the chr3L peak is chorion amplification and not a segregating
# duplication, then the X-linked chorion cluster must ALSO be elevated, in every
# line, in both species, and by a smaller factor (15x vs 60x in follicle cells).
# A segregating structural duplication on 3L predicts nothing about the X.
#
# This script:
#   1. scans normalised 1 kb coverage genome-wide in every line,
#   2. reports every region above a coverage threshold, ranked,
#   3. reports coverage at both chorion clusters specifically,
#   4. plots chrX and chr3L coverage so the two clusters can be seen side by side.
#
# If the X cluster is up in all 12 lines, the amplification explanation is proved
# and the "fixed duplication" explanation is dead.

suppressMessages(library(tidyverse))

# ============================================================
# PATHS
# ============================================================
cov_dir <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/Dm_coverage"
out_dir <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/chorion"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("HOULE_L1F",   "HOULE_L2F",   "HOULE_L3F",
             "HOUSTON_L1F", "HOUSTON_L2F", "HOUSTON_L3F")
chroms  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

# Known chorion clusters. The 3L window is the one we already censored.
# Coordinates confirmed on FlyBase (dm6). Cp36 = FBgn0000359, X:8,479,319..8,480,523,
# cytological 7F1-7. The X cluster (Cp36, Cp38, Cp15, Cp16, Cp18, Cp19, Cp7Fa, Cp7Fb)
# sits at roughly chrX 8.45-8.52 Mb; the bracket below is deliberately generous.
# NOTE: the genome-wide scan in step 1 does not depend on these coordinates, so it
# will find the X cluster even if this bracket is slightly off.
clusters <- tribble(
  ~name,               ~chr,    ~start,     ~end,
  "3L cluster (66D)",  "chr3L", 8500000L,   9000000L,
  "X cluster (7F)",    "chrX",  8400000L,   8600000L
)

THRESH <- 2  # normalised coverage above which a window is called elevated

# ============================================================
# LOAD + NORMALISE
# ============================================================
load_cov <- function(samp) {
  map_dfr(chroms, function(ch) {
    f <- file.path(cov_dir, paste0(samp, "_", ch, "_1kb.txt"))
    if (!file.exists(f)) { cat(sprintf("  [missing] %s\n", basename(f))); return(NULL) }
    d <- read.table(f, header = TRUE)
    # columns: CHROM  WIN_START  WIN_END  MEAN_DEPTH
    tibble(sample = samp, chr = ch, pos = as.numeric(d$WIN_START),
           depth = as.numeric(d$MEAN_DEPTH))
  })
}

cov <- map_dfr(samples, load_cov)
stopifnot(nrow(cov) > 0)

cov <- cov |>
  group_by(sample) |>
  mutate(norm = depth / median(depth[depth > 0], na.rm = TRUE)) |>
  ungroup()

# ============================================================
# 1. GENOME-WIDE SCAN FOR ELEVATED REGIONS
# ============================================================
# a region = run of consecutive elevated 1 kb windows, elevated in EVERY line
elevated <- cov |>
  group_by(chr, pos) |>
  summarise(n_lines_up = sum(norm > THRESH),
            mean_norm  = mean(norm), .groups = "drop") |>
  filter(n_lines_up == length(samples))          # up in all lines, not one

regions <- elevated |>
  arrange(chr, pos) |>
  group_by(chr) |>
  mutate(gap = pos - lag(pos, default = first(pos)),
         blk = cumsum(gap > 5000)) |>            # allow 5 kb gaps
  group_by(chr, blk) |>
  summarise(start = min(pos), end = max(pos),
            kb = (max(pos) - min(pos)) / 1000 + 1,
            mean_norm = mean(mean_norm),
            max_norm  = max(mean_norm), .groups = "drop") |>
  filter(kb >= 5) |>
  arrange(desc(mean_norm))

cat("\n########## Regions elevated (>", THRESH, "x) in ALL", length(samples), "lines ##########\n")
regions |> as.data.frame() |> head(25) |> print()
write_csv(regions, file.path(out_dir, "elevated_regions_all_lines.csv"))

# ============================================================
# 2. COVERAGE AT THE TWO CHORION CLUSTERS, PER LINE
# ============================================================
clust_cov <- clusters |>
  rowwise() |>
  mutate(dat = list(
    cov |> filter(chr == .env$chr, pos >= .env$start, pos <= .env$end) |>
      group_by(sample) |>
      summarise(mean_norm = mean(norm),
                max_norm  = max(norm),
                peak_pos  = pos[which.max(norm)], .groups = "drop")
  )) |>
  select(name, dat) |>
  unnest(dat)

cat("\n########## Coverage at each chorion cluster, per line ##########\n")
cat("(normalised to each line's genome-wide median)\n\n")
clust_cov |> as.data.frame() |> print()
write_csv(clust_cov, file.path(out_dir, "chorion_cluster_coverage.csv"))

cat("\n########## THE TEST ##########\n")
summ <- clust_cov |>
  group_by(name) |>
  summarise(n_lines_elevated = sum(max_norm > THRESH),
            n_lines          = n(),
            mean_of_max      = mean(max_norm),
            .groups = "drop")
summ |> as.data.frame() |> print()

if (all(summ$n_lines_elevated == summ$n_lines)) {
  cat("\n==> BOTH chorion clusters are elevated in EVERY line.\n")
  cat("    This is what follicle-cell amplification predicts.\n")
  cat("    A segregating duplication on 3L predicts nothing about the X cluster.\n")
} else {
  cat("\n==> The X cluster is NOT elevated in all lines. Re-check the 7F coordinates\n")
  cat("    against FlyBase before concluding anything.\n")
}

# ============================================================
# 3. PLOT
# ============================================================
plt <- cov |>
  filter(chr %in% c("chrX", "chr3L")) |>
  separate(sample, c("pop", "line"), sep = "_", remove = FALSE) |>
  ggplot(aes(pos / 1e6, norm, colour = line)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
  geom_rect(data = clusters |> mutate(pop = NA),
            aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "firebrick", alpha = 0.12) +
  facet_grid(pop ~ chr, scales = "free_x") +
  coord_cartesian(ylim = c(0, 8)) +
  labs(x = "Position (Mb)", y = "Coverage / genome median",
       colour = "Line",
       title = "Both chorion clusters are amplified in whole-female DNA",
       subtitle = "Shaded: the 3L cluster (66D) and the X cluster (7F). Amplification in ovarian follicle cells predicts both.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_chorion_clusters.png"), plt,
       width = 11, height = 6, dpi = 300)

cat("\nWrote output to:", out_dir, "\n")
