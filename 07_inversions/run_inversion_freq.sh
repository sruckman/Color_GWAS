#!/bin/bash

#SBATCH --job-name=inversion_freq
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=8G
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

REFALT_DIR="process/Dm_color"
REF="ref/dm6.fa"
OUT_DIR="process/inversions"
SCRIPTS_DIR="scripts"

DROSEU_DIR="tools/DrosEU_pipeline"
MARKER_FILE="process/inversions/inversion_markers_v6.txt"
MARKER_POS="process/inversions/inversion_markers_v6.txt_pos"

SYNC_FILE="${OUT_DIR}/all_populations.sync"
POPS_FILE="${OUT_DIR}/populations.txt"
FILTERED_SYNC="${OUT_DIR}/inversion_markers.sync"
FREQ_FILE="${OUT_DIR}/inversion_frequencies.txt"

mkdir -p "${OUT_DIR}" tools

# ============================================================
# STEP 1 — Clone DrosEU pipeline (if not already present)
# ============================================================

if [ ! -d "${DROSEU_DIR}" ]; then
    echo "=== Cloning DrosEU pipeline ==="
    git clone https://github.com/capoony/DrosEU_pipeline.git "${DROSEU_DIR}"
else
    echo "=== DrosEU pipeline already present at ${DROSEU_DIR} ==="
fi

# ============================================================
# STEP 2 — Convert RefAlt to sync format
# ============================================================

echo ""
echo "=== Converting RefAlt tables to sync format ==="
python3 "${SCRIPTS_DIR}/refalt_to_sync.py" \
    --markers    "${MARKER_FILE}" \
    --refalt-dir "${REFALT_DIR}" \
    --ref        "${REF}" \
    --output     "${SYNC_FILE}" \
    --pops       "${POPS_FILE}"

if [ $? -ne 0 ]; then
    echo "ERROR: refalt_to_sync.py failed"
    exit 1
fi

echo "Sync file: ${SYNC_FILE}"
echo "Populations file: ${POPS_FILE}"
echo ""
cat "${POPS_FILE}"

# ============================================================
# STEP 3 — Sync file is already filtered to marker positions
#           (OverlapSNPs.py skipped — it is Python 2)
# ============================================================

echo ""
echo "=== Sync file already contains only marker positions ==="
echo "  Positions in sync: $(wc -l < ${SYNC_FILE})"

# ============================================================
# STEP 4 — Calculate inversion frequencies (Python 3 replacement for InvFreq.py)
# ============================================================

echo ""
echo "=== Calculating inversion frequencies ==="
python3 "${SCRIPTS_DIR}/calc_inv_freq.py" \
    --sync "${SYNC_FILE}" \
    --inv  "${MARKER_FILE}" \
    --pop  "${POPS_FILE}" \
    > "${FREQ_FILE}"

echo ""
echo "=== Inversion frequencies ==="
cat "${FREQ_FILE}"

echo ""
echo "=== Done ==="
echo "Results: ${FREQ_FILE}"
