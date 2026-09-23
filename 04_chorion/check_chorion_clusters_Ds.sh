#!/bin/bash
#SBATCH --job-name=chorion_Ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%j.out

# D. simulans chorion check. Conserved-cluster replicate of the D. melanogaster proof.
#
# STEP 0 (only if D. simulans 1 kb coverage does not already exist):
#   The extractor already exists at Color/Coverage/get_coverage_Ds.sh.
#   From the Dsim project root, run it against the Dsim BAMs, e.g.:
#     cd /dfs7/adl/sruckman/XQTL_sim/XQTL2
#     sbatch scripts/get_coverage_Ds.sh data/bam/color process/Ds_coverage
#   It is a 1-5 array (one chromosome per task) and writes
#   <sample>_<chr>_1kb.txt with a CHROM/WIN_START/WIN_END/MEAN_DEPTH header,
#   using the downsampled OR original BAMs (either is fine; coverage is normalised
#   per line). Use the reheadered chr-named BAMs if you have them; get_coverage_Ds.sh
#   already handles the NCBI -> chr renaming if you point it at the original BAMs.
#
# STEP 1: this job. First check the coverage files exist and their path matches
#   cov_dir in check_chorion_clusters_Ds.R:
#     ls /dfs7/adl/sruckman/XQTL_sim/XQTL2/process/Ds_coverage/Dsim_1-L1F_chrX_1kb.txt
#   If they are somewhere else, edit cov_dir in the R script.
#
# Upload:
#   scp check_chorion_clusters_Ds.R  sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
#   scp check_chorion_clusters_Ds.sh sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
# Submit:
#   cd /dfs7/adl/sruckman/XQTL/XQTL2 && sbatch scripts/check_chorion_clusters_Ds.sh
# Download:
#   scp -r sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/chorion .
#
# Verdict is printed under "THE TEST (D. simulans)".
# The genome-wide scan is coordinate-free, so it finds the clusters even if the
# chrX bracket in the R script is slightly off. If the X window misses, read the
# elevated chrX region's coordinates from the scan and update the R script.

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/check_chorion_clusters_Ds.R
