#!/bin/bash
#SBATCH --job-name=chorion_check
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=2:00:00
#SBATCH --output=slurm-%j.out

# Proves the chr3L peak is chorion amplification in female DNA, not a duplication.
#
# Drosophila has TWO chorion clusters, both amplified in ovarian follicle cells:
#   chr3L 66D  (~60-fold in follicle cells)  <- the peak we already see
#   chrX  7F   (~15-fold in follicle cells)  <- must ALSO be elevated if this is amplification
#
# A segregating duplication on 3L predicts nothing about the X cluster.
# So if the X cluster is elevated in all 12 lines, the amplification explanation is proved.
#
# Reads the existing 1 kb coverage files in process/Dm_coverage/. No new sequencing.
# Short job; 2 cores x 6G is plenty.
#
# Upload:
#   scp check_chorion_clusters.R  sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
#   scp check_chorion_clusters.sh sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/scripts/
# Submit:
#   cd /dfs7/adl/sruckman/XQTL/XQTL2 && sbatch scripts/check_chorion_clusters.sh
# Download:
#   scp -r sruckman@hpc3.rcic.uci.edu:/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds/plots_censored/chorion .
#
# The verdict is printed to slurm-<jobid>.out under "THE TEST".

module load R/4.2.2

cd /dfs7/adl/sruckman/XQTL/XQTL2

Rscript scripts/check_chorion_clusters.R
