#!/bin/bash
#SBATCH --job-name=null_qq_sim
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1-12:00:00
#SBATCH --output=slurm-%j.out

# simulate_null_qq.R models LD blocks with exponential decay and uses Fisher's
# exact test (USE_FISHER <- TRUE). Keep n_snps = 1e6; do NOT reduce it.
#
# Pass s_mean as the first argument; output filenames are tagged with it, so two
# values can run at once without overwriting each other:
#   sbatch scripts/run_null_qq.sh 0.08
#   sbatch scripts/run_null_qq.sh 0.10
# If no argument is given, the script default (0.08) is used.
#
# Resources: the R script is single-threaded. A full 1e6-SNP / 50-rep Fisher run
# took ~18.75 h and ~4.34 GB. We request 1 core x 5G and 1.5 days wall time.

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/simulate_null_qq.R "$1"
