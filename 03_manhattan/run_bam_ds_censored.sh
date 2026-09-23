#!/bin/bash

#SBATCH --job-name=bam_ds_censored
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00

module load R

Rscript /dfs7/adl/sruckman/XQTL/XQTL2/scripts/plot_bam_ds_censored.R
