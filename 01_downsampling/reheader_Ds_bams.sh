#!/bin/bash

#SBATCH --job-name=reheader_Ds
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00

module load samtools

cd /dfs7/adl/sruckman/XQTL_sim/XQTL2

PIPELINE_SCRIPT="scripts/bam2bcf2REFALT.sh"

# NC -> chr name mapping for D. simulans
# NC_052520.2 = chr2L, NC_052521.2 = chr2R, NC_052522.2 = chr3L
# NC_052523.2 = chr3R, NC_052525.2 = chrX
SED_CMD='s/SN:NC_052520\.2/SN:chr2L/g;
         s/SN:NC_052521\.2/SN:chr2R/g;
         s/SN:NC_052522\.2/SN:chr3L/g;
         s/SN:NC_052523\.2/SN:chr3R/g;
         s/SN:NC_052525\.2/SN:chrX/g'

for comp in Dark_vs_Light Light_vs_Control Dark_vs_Control; do
    echo ""
    echo "============================================"
    echo "Reheadering: ${comp}"
    echo "============================================"

    bam_dir="data/bam/color_ds/${comp}"
    refalt_dir="process/down_sample/bam_ds/${comp}/Ds_color_sep"
    bam_list="helpfiles/color_ds/${comp}/bam_list.txt"

    for bam in "${bam_dir}"/*.bam; do
        [[ "${bam}" == *.bai ]] && continue
        echo "  Reheadering: ${bam}"
        tmp="${bam}.tmp.bam"
        samtools view -H "${bam}" | sed "${SED_CMD}" | \
            samtools reheader - "${bam}" > "${tmp}"
        mv "${tmp}" "${bam}"
        samtools index "${bam}"
        echo "  Done: ${bam}"
    done

    echo "  Re-submitting bam2bcf2REFALT for ${comp}..."
    sbatch --array=1-5 \
        --job-name="refalt_Ds_${comp}" \
        "${PIPELINE_SCRIPT}" \
        "${bam_list}" \
        "${refalt_dir}"
done

echo ""
echo "=== All reheadering done and RefAlt jobs resubmitted ==="
echo "Wait for refalt_Ds_* jobs to finish, then run:"
echo "  sbatch scripts/submit_fisher_Ds_bam_ds.sh"
