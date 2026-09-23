#!/bin/bash
#SBATCH --job-name=chr3L_cnv
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=1:00:00
#SBATCH --output=slurm-%j.out

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/plot_chr3L_cnv.R
