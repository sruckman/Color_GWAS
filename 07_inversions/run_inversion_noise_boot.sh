#!/bin/bash
#SBATCH --job-name=inv_noise_boot
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --time=4:00:00
#SBATCH --output=slurm-%j.out

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/run_inversion_noise_boot.R
