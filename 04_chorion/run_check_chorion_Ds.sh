#!/bin/bash
#SBATCH --job-name=chorion_check_Ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
#SBATCH --output=slurm-%j.out

module load R/4.2.2

Rscript /dfs7/adl/sruckman/XQTL/XQTL2/scripts/check_chorion_clusters_Ds.R
