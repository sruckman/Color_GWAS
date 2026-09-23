#!/bin/bash

#SBATCH --job-name=inv_freq_bam
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00

module load samtools
module load python

export LD_LIBRARY_PATH=/opt/apps/R/4.5.2/lib64/R/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/data/homezvol3/sruckman/.local/lib/python3.9/site-packages:$PYTHONPATH

cd /dfs7/adl/sruckman/XQTL/XQTL2

# ============================================================
# SETTINGS
# ============================================================

BAM_DIR="data/bam/Dm_color"
REF="ref/dm6.fa"
MARKER_FILE="process/inversions/inversion_markers_v6.txt"
BED_FILE="process/inversions/inversion_markers_dm6.bed"
PILEUP_FILE="process/inversions/marker_pileup.txt"
FREQ_FILE="process/inversions/inversion_frequencies_bam.txt"
SCRIPTS_DIR="scripts"

SAMPLES=(HOULE_L1F HOULE_L2F HOULE_L3F HOUSTON_L1F HOUSTON_L2F HOUSTON_L3F)

mkdir -p process/inversions

# ============================================================
# STEP 1 — Run samtools mpileup at all marker positions
# ============================================================

echo "=== Running samtools mpileup at marker positions ==="

# Build BAM list in order
BAMS=()
for s in "${SAMPLES[@]}"; do
    BAMS+=("${BAM_DIR}/${s}.bam")
done

samtools mpileup \
    -l "${BED_FILE}" \
    -f "${REF}" \
    -Q 20 \
    -q 20 \
    -a \
    --no-output-ins \
    --no-output-del \
    --no-output-ends \
    "${BAMS[@]}" \
    > "${PILEUP_FILE}"

echo "  Pileup lines: $(wc -l < ${PILEUP_FILE})"

# ============================================================
# STEP 2 — Parse pileup and calculate inversion frequencies
# ============================================================

echo ""
echo "=== Calculating inversion frequencies from pileup ==="

python3 "${SCRIPTS_DIR}/pileup_to_inv_freq.py" \
    --pileup   "${PILEUP_FILE}" \
    --markers  "${MARKER_FILE}" \
    --samples  "${SAMPLES[@]}" \
    > "${FREQ_FILE}"

echo ""
echo "=== Inversion frequencies (from BAMs) ==="
cat "${FREQ_FILE}"

echo ""
echo "=== Done ==="
echo "Results: ${FREQ_FILE}"
