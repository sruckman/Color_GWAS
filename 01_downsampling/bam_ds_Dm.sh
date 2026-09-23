#!/bin/bash

#SBATCH --job-name=bam_ds_Dm
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00

module load samtools

cd /dfs7/adl/sruckman/XQTL/XQTL2

# ============================================================
# SETTINGS
# ============================================================

BAM_DIR="data/bam/Dm_color"
OUT_BAM_BASE="data/bam/Dm_color_ds"
REFALT_BASE="process/down_sample/bam_ds"
HELPFILES_BASE="helpfiles/color_ds"
PIPELINE_SCRIPT="scripts/bam2bcf2REFALT.sh"
SEED=42

ALL_SAMPLES=(HOULE_L1F HOULE_L2F HOULE_L3F HOUSTON_L1F HOUSTON_L2F HOUSTON_L3F)

declare -A COMP_G1 COMP_G2
COMP_G1[Dark_vs_Light]="HOULE_L3F HOUSTON_L3F"
COMP_G2[Dark_vs_Light]="HOULE_L2F HOUSTON_L2F"
COMP_G1[Light_vs_Control]="HOULE_L2F HOUSTON_L2F"
COMP_G2[Light_vs_Control]="HOULE_L1F HOUSTON_L1F"
COMP_G1[Dark_vs_Control]="HOULE_L3F HOUSTON_L3F"
COMP_G2[Dark_vs_Control]="HOULE_L1F HOUSTON_L1F"

# ============================================================
# STEP 1: COUNT READS FOR ALL SAMPLES
# ============================================================

echo "=== Counting reads for all Dmel samples ==="
declare -A READS
for s in "${ALL_SAMPLES[@]}"; do
    n=$(samtools flagstat "${BAM_DIR}/${s}.bam" | head -1 | awk '{print $1}')
    READS[$s]=$n
    echo "  ${s}: ${n}"
done

# ============================================================
# STEP 2: PER-COMPARISON DOWNSAMPLING + bam2bcf2REFALT
# ============================================================

for comp in Dark_vs_Light Light_vs_Control Dark_vs_Control; do
    echo ""
    echo "============================================"
    echo "Comparison: ${comp}"
    echo "============================================"

    samples=(${COMP_G1[$comp]} ${COMP_G2[$comp]})

    # Find minimum reads among the 4 involved samples
    min_reads=""
    for s in "${samples[@]}"; do
        r=${READS[$s]}
        if [ -z "$min_reads" ] || [ "$r" -lt "$min_reads" ]; then
            min_reads=$r
        fi
    done
    echo "  Target depth (min reads): ${min_reads}"

    out_bam_dir="${OUT_BAM_BASE}/${comp}"
    refalt_dir="${REFALT_BASE}/${comp}/Dm_color_sep"
    helpfiles_dir="${HELPFILES_BASE}/${comp}"
    mkdir -p "${out_bam_dir}" "${refalt_dir}" "${helpfiles_dir}"

    bam_list="${helpfiles_dir}/bam_list.txt"
    > "${bam_list}"

    for s in "${samples[@]}"; do
        r=${READS[$s]}
        in_bam="${BAM_DIR}/${s}.bam"
        out_bam="${out_bam_dir}/${s}.bam"

        if [ "${r}" -le "${min_reads}" ]; then
            echo "  ${s}: no downsampling needed (${r} <= ${min_reads}), copying"
            cp "${in_bam}" "${out_bam}"
        else
            seed_frac=$(python3 -c "
frac = ${min_reads} / ${r}
decimal = int(frac * 1000000)
print(f'${SEED}.{decimal:06d}')
")
            echo "  ${s}: ${r} -> ${min_reads} reads (fraction ${seed_frac})"
            samtools view -s "${seed_frac}" -b -@ 3 "${in_bam}" -o "${out_bam}"
        fi

        samtools index "${out_bam}"
        echo "${out_bam}" >> "${bam_list}"
    done

    echo "  bam_list written: ${bam_list}"
    cat "${bam_list}"

    echo "  Submitting bam2bcf2REFALT for ${comp}..."
    sbatch --array=1-5 \
        --job-name="refalt_Dm_${comp}" \
        "${PIPELINE_SCRIPT}" \
        "${bam_list}" \
        "${refalt_dir}"
done

echo ""
echo "=== All comparisons submitted ==="
echo "Wait for refalt_Dm_* jobs to finish, then run:"
echo "  sbatch scripts/submit_fisher_Dm_bam_ds.sh"
