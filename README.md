# Color_GWAS

Analysis code for:

> Selection on cuticle color produces a genome-wide response that does not localize to known major genes in *Drosophila*. (In preparation.)

We selected on thoracic trident color for 16 generations in two populations of *Drosophila melanogaster* (Toronto and Miami) and two of *D. simulans* (Tallahassee-1 and Tallahassee-2). At generation 16 we pool-sequenced the dark, light, and control lines and tested each SNP for allele-frequency differentiation. This repository holds the code for every analysis and figure in the paper.

## Before you run anything

- **Cluster paths are hardcoded.** Scripts were run on UC Irvine's HPC3 cluster and point to `/dfs7/adl/sruckman/...`. Edit these paths for your own system.
- **Folder layout differs from the cluster.** On the cluster, every script sat in one flat `scripts/` folder, and jobs were submitted from the project root. Here the scripts are grouped by analysis step. Some scripts call others by `scripts/<name>`, so copy them into one `scripts/` folder or edit those calls.
- **SLURM settings are lab-specific.** The `.sh` wrappers use `#SBATCH -A tdlong_lab -p standard`. Change the account and partition for your cluster.

## Population codes

The scripts use internal codes. The paper uses place names.

| Script code | Paper name | Species |
|---|---|---|
| `HOULE` | Toronto | *D. melanogaster* |
| `HOUSTON` | Miami | *D. melanogaster* |
| `Dsim_1` | Tallahassee-1 | *D. simulans* |
| `Dsim_2` | Tallahassee-2 | *D. simulans* |

Lines: `L1F` = control, `L2F` = light selected, `L3F` = dark selected.

## Upstream pipeline

Read alignment and allele counting used the Long lab XQTL2 pipeline: https://github.com/tdlong/XQTL2. That covers BWA-MEM alignment, read groups, and `bam2bcf2REFALT.sh`, which turns BAMs into per-chromosome RefAlt allele-count tables. Everything in this repository starts from those BAMs and RefAlt files.

## Run order

| Step | Folder | What it does |
|---|---|---|
| 1 | `01_downsampling/` | Downsamples BAMs within each contrast to equal coverage (`samtools view -s`, seed 42). Renames *D. simulans* NCBI chromosome names to chr-style names. Submits the XQTL2 RefAlt jobs. |
| 2 | `02_fisher/` | Per-SNP Fisher's exact tests of allele-frequency differentiation for each contrast. `count_snps.sh` counts the SNPs tested per contrast. |
| 3 | `03_manhattan/` | QQ and Manhattan plots, with and without the chromosome 3L chorion region censored. |
| 4 | `04_chorion/` | 1 kb coverage along each chromosome and coverage over the two chorion gene clusters (66D on 3L, 7F on X). Tests whether the 3L peak comes from follicle-cell chorion amplification. |
| 5 | `05_null_simulation/` | Wright-Fisher forward simulation of drift only and of variable polygenic selection. Builds the null QQ envelopes. |
| 6 | `06_power_model/` | Power model comparing our E&R design (R = 1, 2, 5, 10 replicate populations) against a single-generation case/control design. |
| 7 | `07_inversions/` | Cosmopolitan inversion frequencies from diagnostic marker SNPs, and tests of whether inversions inflate background noise (block bootstrap). |
| 8 | `08_fst/` | Hudson's Fst between the two control lines within each species, in 10 kb and 50 kb windows. |

## Figures

| Figure | Script |
|---|---|
| Figure 1 (QQ plots with null envelopes) | `05_null_simulation/simulate_null_qq.R` |
| Figures 2 and 3 (Manhattan plots) | `03_manhattan/plot_bam_ds_manhattan.R` (censored versions: `plot_bam_ds_censored.R`) |
| Figures 4 and 5 (chorion amplification) | `04_chorion/plot_chorion_evidence.R` |
| Figure 6 (power) | `06_power_model/make_figures_Dm.R`, `make_figures_Ds.R` |
| Figure 7 (inversion frequencies) | `07_inversions/plot_inversion_freq.R` |
| Figure 8 (Fst) | `08_fst/run_fst_Dm_Ds.R` |
| Figure S1 (noise by starting frequency) | `06_power_model/make_figures_Dm.R`, `make_figures_Ds.R` |

## Folder notes

**`06_power_model/`.** `design_power_model.R` holds the core power-model functions and all parameters, written by A. D. Long (Long et al. 2026). `run_design_power_Dm.R` and `run_design_power_Ds.R` run the model on each species' data. `make_figures_Dm.R` and `make_figures_Ds.R` source `design_power_model.R` and make the power and noise figures with block-bootstrap intervals.

**`07_inversions/`.**
- `inversion_markers_dm3.bed` lists the diagnostic marker SNPs in dm3 coordinates. `liftover_markers.sh` and `remap_markers.py` move them to dm6.
- `run_inv_freq_bam.sh` and `pileup_to_inv_freq.py` estimate inversion frequencies from BAM pileups. This is the method used in the paper.
- `run_inversion_freq.sh`, `refalt_to_sync.py`, and `calc_inv_freq.py` are an earlier RefAlt-based route to the same estimates.
- `inversion_frequencies_bam.txt` and `marker_pileup_counts.txt` are the output that `plot_inversion_freq.R` reads.
- `run_inversion_noise_Dm.R` and `run_inversion_noise_boot.R` run the noise comparison and its bootstrap.

**Censored regions.** Some analyses exclude the chorion peak on chromosome 3L. The censored ranges are set at the top of each script.

## Software

- R 4.2.2 with tidyverse
- Python 3 (HPC3 `python` module)
- samtools and BWA (HPC3 modules)
- SLURM

## License

MIT. See `LICENSE`.
