#!/bin/bash

#SBATCH --job-name=fisher_Ds_bam_ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=2
#SBATCH --time=8:00:00
#SBATCH --array=1-5

module load R

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

echo "Starting BAM-downsampled Dsim GWAS for chromosome: $mychr"
Rscript /dfs7/adl/sruckman/XQTL_sim/XQTL2/scripts/run_fisher_Ds_bam_ds.R "$mychr"
echo "Finished chromosome: $mychr"
