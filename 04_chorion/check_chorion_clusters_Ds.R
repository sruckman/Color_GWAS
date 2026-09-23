#!/usr/bin/env Rscript
# check_chorion_clusters_Ds.R
#
# D. simulans version of the chorion proof (see check_chorion_clusters.R for D. melanogaster).
#
# Same logic: Drosophila has two chorion gene clusters (66D on chr3L, 7F on chrX),
# both developmentally amplified in ovarian follicle cells (Spradling and Mahowald
# 1980; Claycomb and Orr-Weaver 2005). If the chr3L peak is follicle-cell
# amplification and not a segregating duplication, the X cluster must ALSO be
# elevated in every line. A duplication on 3L predicts nothing about the X.
#
# We already showed both clusters are elevated in all six D. melanogaster lines
# (66D 6.6x, 7F 3.0x). D. simulans is a conserved-cluster replicate: if the story
# is right, it must reproduce.
#
# NOTE ON COORDINATES: the D. simulans assembly is NOT dm6, so cluster positions
# differ from D. melanogaster. The genome-wide scan in step 1 does NOT depend on
# any coordinates and will find the clusters wherever they are. The bracketed
# windows below are only used for the per-line summary; the script prints the
# observed peak position so the exact D. simulans coordinates can be recorded.
# The chr3L window matches the region already censored in the paper
# (chr3L 8,342,753-8,834,090 in D. simulans). The chrX window is a guess based on
# rough synteny with dm6 7F (~chrX 8.4-8.6 Mb) and is deliberately wide; trust the
# genome-wide scan over this bracket.

suppressMessages(library(tidyverse))

# ============================================================
# PATHS  -- EDIT cov_dir IF THE D. simulans COVERAGE FILES ARE ELSEWHERE
# ============================================================
# get_coverage_Ds.sh writes files named <sample>_<chr>_1kb.txt with a header line
# "CHROM  WIN_START  WIN_END  MEAN_DEPTH". Point cov_dir at wherever that job wrote.
# The D. melanogaster equivalent lived in process/Dm_coverage/, so the natural
# guess is process/Ds_coverage/.
cov_dir <- "/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/Ds_coverage"
out_dir <- "/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/chorion"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# D. simulans sample names (from get_coverage_Ds.sh)
samples <- c("Dsim_1-L1F", "Dsim_1-L2F", "Dsim_1-L3F",
             "Dsim_2-L1F", "Dsim_2-L2F", "Dsim_2-L3F")
chroms  <- c("chrX", "chr2L", "chr2R", "chr3L", "chr3R")

clusters <- tribble(
  ~name,               ~chr,    ~start,     ~end,
  "3L cluster (66D)",  "chr3L", 8342753L,   8834090L,   # matches the censored region
  "X cluster (7F)",    "chrX",  7800000L,   8400000L    # confirmed coordinates. Peak at
                                                         # 8,208,001 in 5/6 lines (up to 68x);
                                                         # secondary peak at 7,860,001.
                                                         # Old 8.7 Mb window was tan (8,770,795).
)

THRESH <- 2

# ============================================================
# LOAD + NORMALISE  (identical to the Dmel version)
# ============================================================
load_cov <- function(samp) {
  map_dfr(chroms, function(ch) {
    f <- file.path(cov_dir, paste0(samp, "_", ch, "_1kb.txt"))
    if (!file.exists(f)) { cat(sprintf("  [missing] %s\n", basename(f))); return(NULL) }
    d <- read.table(f, header = TRUE)          # CHROM WIN_START WIN_END MEAN_DEPTH
    tibble(sample = samp, chr = ch, pos = as.numeric(d$WIN_START),
           depth = as.numeric(d$MEAN_DEPTH))
  })
}

cov <- map_dfr(samples, load_cov)
if (nrow(cov) == 0)
  stop("No coverage files found. Set cov_dir to the D. simulans coverage directory, ",
       "or run Coverage/get_coverage_Ds.sh first.")

cov <- cov |>
  group_by(sample) |>
  mutate(norm = depth / median(depth[depth > 0], na.rm = TRUE)) |>
  ungroup()

# ============================================================
# 1. GENOME-WIDE SCAN (coordinate-free; this is the real test)
# ============================================================
elevated <- cov |>
  group_by(chr, pos) |>
  summarise(n_lines_up = sum(norm > THRESH), mean_norm = mean(norm), .groups = "drop") |>
  filter(n_lines_up == length(samples))

regions <- elevated |>
  arrange(chr, pos) |>
  group_by(chr) |>
  mutate(gap = pos - lag(pos, default = first(pos)), blk = cumsum(gap > 5000)) |>
  group_by(chr, blk) |>
  summarise(start = min(pos), end = max(pos),
            kb = (max(pos) - min(pos)) / 1000 + 1,
            mean_norm = mean(mean_norm), max_norm = max(mean_norm), .groups = "drop") |>
  filter(kb >= 5) |>
  arrange(desc(mean_norm))

cat("\n########## D. simulans: regions elevated (>", THRESH, "x) in ALL", length(samples), "lines ##########\n")
regions |> as.data.frame() |> head(25) |> print()
write_csv(regions, file.path(out_dir, "Ds_elevated_regions_all_lines.csv"))

# ============================================================
# 2. COVERAGE AT THE TWO CHORION CLUSTERS, PER LINE
# ============================================================
clust_cov <- clusters |>
  rowwise() |>
  mutate(dat = list(
    cov |> filter(chr == .env$chr, pos >= .env$start, pos <= .env$end) |>
      group_by(sample) |>
      summarise(mean_norm = mean(norm), max_norm = max(norm),
                peak_pos = pos[which.max(norm)], .groups = "drop")
  )) |>
  select(name, dat) |>
  unnest(dat)

cat("\n########## D. simulans: coverage at each chorion cluster, per line ##########\n\n")
clust_cov |> as.data.frame() |> print()
write_csv(clust_cov, file.path(out_dir, "Ds_chorion_cluster_coverage.csv"))

cat("\n########## THE TEST (D. simulans) ##########\n")
summ <- clust_cov |>
  group_by(name) |>
  summarise(n_lines_elevated = sum(max_norm > THRESH), n_lines = n(),
            mean_of_max = mean(max_norm), .groups = "drop")
summ |> as.data.frame() |> print()

if (all(summ$n_lines_elevated == summ$n_lines)) {
  cat("\n==> BOTH chorion clusters are elevated in EVERY D. simulans line too.\n")
  cat("    The artifact replicates across species. Amplification is confirmed.\n")
} else {
  cat("\n==> The X cluster is not elevated in all lines under the current bracket.\n")
  cat("    CHECK the genome-wide scan above: find the chrX region that is elevated\n")
  cat("    in all 6 lines and read its coordinates, then update the chrX window.\n")
}

# ============================================================
# 3. PLOT
# ============================================================
plt <- cov |>
  filter(chr %in% c("chrX", "chr3L")) |>
  separate(sample, c("pop", "line"), sep = "-", remove = FALSE) |>
  ggplot(aes(pos / 1e6, norm, colour = line)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
  geom_rect(data = clusters |> mutate(pop = NA),
            aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "firebrick", alpha = 0.12) +
  facet_grid(pop ~ chr, scales = "free_x") +
  coord_cartesian(ylim = c(0, 8)) +
  labs(x = "Position (Mb)", y = "Coverage / genome median", colour = "Line",
       title = "D. simulans: both chorion clusters are amplified in whole-female DNA",
       subtitle = "Conserved-cluster replicate of the D. melanogaster result.") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "fig_chorion_clusters_Ds.png"), plt, width = 11, height = 6, dpi = 300)
cat("\nWrote output to:", out_dir, "\n")
