#!/bin/bash
# Liftover inversion marker SNPs from dm3 to dm6
# Run from /dfs7/adl/sruckman/XQTL/XQTL2

set -e

TOOLS_DIR="tools"
BED_IN="scripts/inversion_markers_dm3.bed"
BED_OUT="process/inversions/inversion_markers_dm6.bed"
UNMAPPED="process/inversions/inversion_markers_dm3_unmapped.bed"
CHAIN="${TOOLS_DIR}/dm3ToDm6.over.chain.gz"
LIFTOVER="${TOOLS_DIR}/liftOver"
FINAL_OUT="process/inversions/inversion_markers_v6.txt"
FINAL_POS="process/inversions/inversion_markers_v6.txt_pos"

mkdir -p process/inversions "${TOOLS_DIR}"

# Download liftOver binary if not present
if [ ! -f "${LIFTOVER}" ]; then
    echo "Downloading liftOver..."
    wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver \
        -O "${LIFTOVER}"
    chmod +x "${LIFTOVER}"
fi

# Download chain file if not present
if [ ! -f "${CHAIN}" ]; then
    echo "Downloading dm3->dm6 chain file..."
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz \
        -O "${CHAIN}"
fi

# Run liftOver
echo "Running liftOver..."
"${LIFTOVER}" "${BED_IN}" "${CHAIN}" "${BED_OUT}" "${UNMAPPED}"

TOTAL=$(wc -l < "${BED_IN}")
LIFTED=$(wc -l < "${BED_OUT}")
LOST=$(grep -v "^#" "${UNMAPPED}" | wc -l)
echo "  Total markers:   ${TOTAL}"
echo "  Lifted to dm6:   ${LIFTED}"
echo "  Unmapped:        ${LOST}"

# Convert BED output to inversion_markers_v6.txt format
# BED name field is: Inversion|chr|Allele
# Output: Inversion <tab> chr (without chr prefix) <tab> pos (1-based) <tab> Allele
echo "Writing inversion_markers_v6.txt..."
awk '{
    split($4, a, "|")
    inv    = a[1]
    allele = a[3]
    chr    = $1
    sub("chr", "", chr)      # strip chr prefix -> 2L, 2R etc.
    pos    = $3              # BED end = 1-based position
    print inv "\t" chr "\t" pos "\t" allele
}' "${BED_OUT}" | sort -k2,2 -k3,3n > "${FINAL_OUT}"

# Write position-only file for OverlapSNPs.py
awk '{print $2 "\t" $3}' "${FINAL_OUT}" > "${FINAL_POS}"

echo ""
echo "=== Done ==="
echo "Marker file (v6): ${FINAL_OUT}  ($(wc -l < ${FINAL_OUT}) lines)"
echo "Position file:    ${FINAL_POS}"
echo ""
echo "Check a few lines:"
head -5 "${FINAL_OUT}"
