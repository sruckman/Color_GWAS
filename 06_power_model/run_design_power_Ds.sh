#!/bin/bash
#SBATCH --job-name=design_power_Ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=3:00:00
#SBATCH --output=slurm-%j.out

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/run_design_power_Ds.R
