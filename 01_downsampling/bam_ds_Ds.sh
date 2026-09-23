#!/bin/bash

#SBATCH --job-name=bam_ds_Ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --time=2-00:00:00

module load samtools

cd /dfs7/adl/sruckman/XQTL_sim/XQTL2

# ============================================================
# SETTINGS
# ============================================================

BAM_DIR="data/bam/color"
OUT_BAM_BASE="data/bam/color_ds"
REFALT_BASE="process/down_sample/bam_ds"
HELPFILES_BASE="helpfiles/color_ds"
PIPELINE_SCRIPT="scripts/bam2bcf2REFALT.sh"
SEED=42

# NC -> chr name mapping for D. simulans BAM headers
SED_CMD='s/SN:NC_052520\.2/SN:chr2L/g;
         s/SN:NC_052521\.2/SN:chr2R/g;
         s/SN:NC_052522\.2/SN:chr3L/g;
         s/SN:NC_052523\.2/SN:chr3R/g;
         s/SN:NC_052525\.2/SN:chrX/g'

ALL_SAMPLES=("Dsim_1-L1F" "Dsim_1-L2F" "Dsim_1-L3F" "Dsim_2-L1F" "Dsim_2-L2F" "Dsim_2-L3F")

declare -A COMP_G1 COMP_G2
COMP_G1[Dark_vs_Light]="Dsim_1-L3F Dsim_2-L3F"
COMP_G2[Dark_vs_Light]="Dsim_1-L2F Dsim_2-L2F"
COMP_G1[Light_vs_Control]="Dsim_1-L2F Dsim_2-L2F"
COMP_G2[Light_vs_Control]="Dsim_1-L1F Dsim_2-L1F"
COMP_G1[Dark_vs_Control]="Dsim_1-L3F Dsim_2-L3F"
COMP_G2[Dark_vs_Control]="Dsim_1-L1F Dsim_2-L1F"

# ============================================================
# STEP 0: CHECK BAM HEADERS (chromosome naming)
# ============================================================

echo "=== Checking Dsim BAM header (chromosome names) ==="
samtools view -H "${BAM_DIR}/Dsim_1-L1F.bam" | grep "^@SQ" | awk '{print $2, $3}' | head -10
echo "(Confirm these match the chromosome names expected by bam2bcf2REFALT)"
echo ""

# ============================================================
# STEP 1: COUNT READS FOR ALL SAMPLES
# ============================================================

echo "=== Counting reads for all Dsim samples ==="
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

    min_reads=""
    for s in "${samples[@]}"; do
        r=${READS[$s]}
        if [ -z "$min_reads" ] || [ "$r" -lt "$min_reads" ]; then
            min_reads=$r
        fi
    done
    echo "  Target depth (min reads): ${min_reads}"

    out_bam_dir="${OUT_BAM_BASE}/${comp}"
    refalt_dir="${REFALT_BASE}/${comp}/Ds_color_sep"
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

        # Reheader: replace NCBI chromosome names with chr-style names
        tmp="${out_bam}.tmp.bam"
        samtools view -H "${out_bam}" | sed "${SED_CMD}" | \
            samtools reheader - "${out_bam}" > "${tmp}"
        mv "${tmp}" "${out_bam}"

        samtools index "${out_bam}"
        echo "${out_bam}" >> "${bam_list}"
    done

    echo "  bam_list written: ${bam_list}"
    cat "${bam_list}"

    echo "  Submitting bam2bcf2REFALT for ${comp}..."
    sbatch --array=1-5 \
        --job-name="refalt_Ds_${comp}" \
        "${PIPELINE_SCRIPT}" \
        "${bam_list}" \
        "${refalt_dir}"
done

echo ""
echo "=== All comparisons submitted ==="
echo "Wait for refalt_Ds_* jobs to finish, then run:"
echo "  sbatch scripts/submit_fisher_Ds_bam_ds.sh"
