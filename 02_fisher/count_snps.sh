#!/bin/bash
# count_snps.sh
# Count SNPs tested per population-by-contrast combination.
# Run interactively on the cluster (no SLURM needed — finishes in seconds).
#
# Usage:
#   bash scripts/count_snps.sh
#
# POS is assumed to be column 2 in the space-delimited GWAS_results files.

DM="/dfs7/adl/sruckman/XQTL/XQTL2/process/down_sample/bam_ds"
DS="/dfs7/adl/sruckman/XQTL_sim/XQTL2/process/down_sample/bam_ds"
CHROMS="chrX chr2L chr2R chr3L chr3R"
COMPS="Dark_vs_Control Light_vs_Control"

# Chorion exclusion windows (per species, per chr).
# Returns awk condition to KEEP rows (i.e., NOT in the chorion window).
chorion_awk_dmel() {
    local chr=$1
    case $chr in
        chr3L) echo '$2 < 8500000 || $2 > 9000000' ;;
        chrX)  echo '$2 < 8400000 || $2 > 8600000' ;;
        *)     echo '1' ;;
    esac
}
chorion_awk_dsim() {
    local chr=$1
    case $chr in
        chr3L) echo '$2 < 8342753 || $2 > 8834090' ;;
        chrX)  echo '$2 < 7800000 || $2 > 8400000' ;;
        *)     echo '1' ;;
    esac
}

printf "\n%-8s %-12s %-22s %12s %16s\n" \
    "Species" "Pop" "Contrast" "Total_SNPs" "Filtered_SNPs"
printf "%-8s %-12s %-22s %12s %16s\n" \
    "-------" "-----------" "---------------------" "----------" "--------------"

for comp in $COMPS; do
    # D. melanogaster
    for pop in HOULE HOUSTON; do
        total=0; filtered=0
        for chr in $CHROMS; do
            f="${DM}/${comp}/Dm_color_sep/GWAS_results_${pop}_${comp}_${chr}.txt"
            if [ ! -f "$f" ]; then echo "MISSING: $f" >&2; continue; fi
            t=$(tail -n +2 "$f" | wc -l)
            cond=$(chorion_awk_dmel $chr)
            filt=$(tail -n +2 "$f" | awk "$cond" | wc -l)
            total=$((total + t))
            filtered=$((filtered + filt))
        done
        printf "%-8s %-12s %-22s %12s %16s\n" \
            "Dmel" "$pop" "$comp" "$total" "$filtered"
    done

    # D. simulans
    for pop in Dsim_1 Dsim_2; do
        total=0; filtered=0
        for chr in $CHROMS; do
            f="${DS}/${comp}/Ds_color_sep/GWAS_results_${pop}_${comp}_${chr}.txt"
            if [ ! -f "$f" ]; then echo "MISSING: $f" >&2; continue; fi
            t=$(tail -n +2 "$f" | wc -l)
            cond=$(chorion_awk_dsim $chr)
            filt=$(tail -n +2 "$f" | awk "$cond" | wc -l)
            total=$((total + t))
            filtered=$((filtered + filt))
        done
        printf "%-8s %-12s %-22s %12s %16s\n" \
            "Dsim" "$pop" "$comp" "$total" "$filtered"
    done

    printf "\n"
done
