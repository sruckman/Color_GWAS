#!/bin/bash
#SBATCH --job-name=fst_Dm_Ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%j.out

# Fst between the two control pools per species.
# D. mel: Toronto (HOULE_L1F) vs Miami (HOUSTON_L1F)
# D. sim: Tallahassee-1 (Dsim_1-L1F) vs Tallahassee-2 (Dsim_2-L1F)
#
# Upload:
#   scp Color/Down_sampled/design_power/run_fst_Dm_Ds.R  sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
#   scp Color/Down_sampled/design_power/run_fst_Dm_Ds.sh sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
# Submit:
#   cd /dfs7/adl/sruckman/XQTL/XQTL2 && sbatch scripts/run_fst_Dm_Ds.sh
# Download:
#   scp -r sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/fst .

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/run_fst_Dm_Ds.R
