#!/bin/bash

# --- Extract 1kb binned coverage from BAM files for D. simulans ---
# Usage: sbatch get_coverage_Ds.sh <bam_dir> <output_dir>
# BAMs are expected to be named <sample>.bam (e.g. Dsim_1-L1F.bam)

# --- SLURM Configuration ---
#SBATCH --job-name=coverage_Ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-5

module load samtools

# NCBI chromosome names in the BAM headers, mapped to simple names
declare -a chrs_ncbi=("NC_052520.2" "NC_052521.2" "NC_052522.2" "NC_052523.2" "NC_052525.2")
declare -a chrs_simple=("chr2L" "chr2R" "chr3L" "chr3R" "chrX")

mychr_ncbi=${chrs_ncbi[$SLURM_ARRAY_TASK_ID - 1]}
mychr_simple=${chrs_simple[$SLURM_ARRAY_TASK_ID - 1]}

bam_dir=$1
output_dir=$2

samples=("Dsim_1-L1F" "Dsim_1-L2F" "Dsim_1-L3F" "Dsim_2-L1F" "Dsim_2-L2F" "Dsim_2-L3F")

mkdir -p $output_dir

for sample in "${samples[@]}"; do
    bam="${bam_dir}/${sample}.bam"
    output="${output_dir}/${sample}_${mychr_simple}_1kb.txt"

    echo -e "CHROM\tWIN_START\tWIN_END\tMEAN_DEPTH" > $output

    # Use NCBI name for samtools depth, rename to simple name in output
    samtools depth -a -r $mychr_ncbi $bam | awk -v chr=$mychr_simple '
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

    echo "Saved coverage for $sample on $mychr_simple to $output"
done

echo "Coverage extraction complete for $mychr_simple"
