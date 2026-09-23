#!/bin/bash
#SBATCH --job-name=inv_noise
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=4:00:00
#SBATCH --output=slurm-%j.out

# Tony's "accidental inversion experiment":
# D. mel-2 (HOUSTON) carries In(2L)t, In(3L)P, In(3R)K.  D. mel-1 (HOULE) carries none.
# Asks whether inversions inflate the genome-wide noise in allele-frequency change,
# both between populations and between inverted vs non-inverted arms.
#
# Reads the same downsampled RefAlt files as run_design_power_Dm.R.
# Runtime is short (it only recomputes sigma(p0) on subsets), 2 cores x 6G is plenty.
#
# Upload:
#   scp run_inversion_noise_Dm.R  sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
#   scp run_inversion_noise_Dm.sh sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
# Submit:
#   cd /dfs7/adl/sruckman/XQTL/XQTL2 && sbatch scripts/run_inversion_noise_Dm.sh
# Download:
#   scp -r sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/inversion_noise .
#
# The key numbers are printed to slurm-<jobid>.out.

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/run_inversion_noise_Dm.R
