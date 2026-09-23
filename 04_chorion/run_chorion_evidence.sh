#!/bin/bash
#SBATCH --job-name=chorion_evidence
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%j.out

# 3-panel chorion evidence figure (chr3L and chrX, per species per contrast).
# Reads GWAS results, 1 kb coverage files, and RefAlt files — no new cluster jobs needed.
# Output: process/down_sample/bam_ds/plots_censored/chorion_evidence/
#
# Upload:
#   scp Color/Down_sampled/manhattan/plot_chorion_evidence.R  sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
#   scp Color/Down_sampled/manhattan/run_chorion_evidence.sh  sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
#
# Submit (from Dmel project root):
#   cd /dfs7/adl/sruckman/XQTL/XQTL2 && sbatch scripts/run_chorion_evidence.sh
#
# Download when done:
#   scp -r sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/chorion_evidence .

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/plot_chorion_evidence.R
