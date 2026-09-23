#!/bin/bash

#SBATCH --job-name=remap_markers
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=2:00:00

module load samtools
module load bwa
module load python

export LD_LIBRARY_PATH=/opt/apps/R/4.5.2/lib64/R/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/data/homezvol3/sruckman/.local/lib/python3.9/site-packages:$PYTHONPATH

cd /dfs7/adl/sruckman/XQTL/XQTL2

# ============================================================
# SETTINGS
# ============================================================

REF_DIR="ref"
DM3_FA="${REF_DIR}/dm3.fa"
DM6_FA="${REF_DIR}/dm6.fa"
SCRIPTS_DIR="scripts"

# dm3-coordinate BED (from liftover_markers.sh step — original Kapun 2014 positions)
DM3_BED="process/inversions/inversion_markers_dm3.bed"

# dm3-format tab file that remap_markers.py reads: inv\tchrom\tpos\tallele
DM3_TAB="process/inversions/inversion_markers_dm3.txt"

# liftOver output BED (produced by liftover_markers.sh) — used for comparison only
LIFTOVER_BED="process/inversions/inversion_markers_dm6_liftover.bed"

OUT_DIR="process/inversions"

mkdir -p "${REF_DIR}" "${OUT_DIR}"

# ============================================================
# STEP 1 — Convert dm3 BED to tab format for remap_markers.py
#   BED cols: chr2L  pos-1  pos  Inv|chrom|Allele
#   Tab cols: Inv    chrom  pos  Allele  (chrom = 2L without chr prefix)
# ============================================================

echo "=== Converting dm3 BED to tab format ==="
awk 'BEGIN{OFS="\t"} {
    n = split($4, a, "|")
    # a[1]=Inv  a[2]=chrom(2L)  a[3]=Allele
    # $3 = BED end = 1-based position
    print a[1], a[2], $3, a[3]
}' "${DM3_BED}" > "${DM3_TAB}"
echo "  Markers: $(wc -l < ${DM3_TAB})"

# ============================================================
# STEP 2 — Download and index dm3 (UCSC dm3 = FlyBase r5)
# ============================================================

if [ ! -f "${DM3_FA}" ]; then
    echo "=== Downloading dm3 reference (this may take a few minutes) ==="
    wget -q -O "${DM3_FA}.gz" \
        "http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/dm3.fa.gz"
    echo "  Decompressing dm3..."
    gunzip "${DM3_FA}.gz"
else
    echo "=== dm3 already present at ${DM3_FA} ==="
fi

if [ ! -f "${DM3_FA}.fai" ]; then
    echo "=== Indexing dm3 with samtools faidx ==="
    samtools faidx "${DM3_FA}"
else
    echo "=== dm3 faidx index already present ==="
fi

# ============================================================
# STEP 3 — Check dm6 BWA index (index if missing)
# ============================================================

if [ ! -f "${DM6_FA}.bwt" ]; then
    echo "=== BWA index not found for dm6 — indexing now (~10 min) ==="
    bwa index "${DM6_FA}"
else
    echo "=== dm6 BWA index already present ==="
fi

# ============================================================
# STEP 4 — Run sequence-based remapping
# ============================================================

LIFTOVER_ARG=""
if [ -f "${LIFTOVER_BED}" ]; then
    echo "=== liftOver BED found — will compare coordinates ==="
    LIFTOVER_ARG="--liftover ${LIFTOVER_BED}"
else
    echo "=== No liftOver BED at ${LIFTOVER_BED} — skipping comparison ==="
fi

echo ""
echo "=== Running sequence-based marker remapping ==="
python3 "${SCRIPTS_DIR}/remap_markers.py" \
    --markers  "${DM3_TAB}" \
    --dm3      "${DM3_FA}" \
    --dm6      "${DM6_FA}" \
    ${LIFTOVER_ARG} \
    --out-dir  "${OUT_DIR}"

echo ""
echo "=== Remap report (first 20 lines) ==="
head -20 "${OUT_DIR}/remap_report.txt"

echo ""
echo "=== Summary ==="
echo "Remapped markers : ${OUT_DIR}/inversion_markers_v6_remapped.txt"
echo "Full report      : ${OUT_DIR}/remap_report.txt"
echo "=== Done ==="
