#!/bin/bash

# --- Extract 1kb binned coverage from BAM files for D. melanogaster ---
# Usage: sbatch get_coverage_Dm.sh <bam_dir> <output_dir>
# BAMs are expected to be named <sample>.bam (e.g. HOULE_L1F.bam)

# --- SLURM Configuration ---
#SBATCH --job-name=coverage_Dm
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5

module load samtools

declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}

bam_dir=$1
output_dir=$2

samples=("HOULE_L1F" "HOULE_L2F" "HOULE_L3F" "HOUSTON_L1F" "HOUSTON_L2F" "HOUSTON_L3F")

mkdir -p $output_dir

for sample in "${samples[@]}"; do
    bam="${bam_dir}/${sample}.bam"
    output="${output_dir}/${sample}_${mychr}_1kb.txt"

    echo -e "CHROM\tWIN_START\tWIN_END\tMEAN_DEPTH" > $output

    samtools depth -a -r $mychr $bam | awk -v chr=$mychr '
    {
        bin = int(($2 - 1) / 1000);
        if (NR == 1) {
            prev_bin = bin; sum = $3; count = 1;
        } else if (bin != prev_bin) {
            printf "%s\t%d\t%d\t%.4f\n", chr, prev_bin*1000+1, (prev_bin+1)*1000, sum/count;
            prev_bin = bin; sum = $3; count = 1;
        } else {
            sum += $3; count++;
        }
    }
    END {
        if (count > 0)
            printf "%s\t%d\t%d\t%.4f\n", chr, prev_bin*1000+1, (prev_bin+1)*1000, sum/count;
    }' >> $output

    echo "Saved coverage for $sample on $mychr to $output"
done

echo "Coverage extraction complete for $mychr"
