#!/usr/bin/env python3
"""
refalt_to_sync.py

Convert RefAlt tables to PoPoolation2 sync format, filtered to inversion marker positions.

Sync format (tab-separated):
    chr  pos  ref_base  A:T:C:G:N:del  A:T:C:G:N:del  ...
    (one count column per sample, in the order they appear in the RefAlt header)

Chromosome naming: outputs 2L/2R/3L/3R/X style to match the DrosEU marker file.
RefAlt files use chr2L style; dm6.fa uses chr2L style — conversion is handled internally.

Usage:
    python3 refalt_to_sync.py \\
        --markers  DrosEU_pipeline/data/inversion_markers_v6.txt \\
        --refalt-dir  process/Dm_color \\
        --ref  ref/dm6.fa \\
        --output  process/inversions/all_populations.sync \\
        --pops   process/inversions/populations.txt
"""

import argparse
import os
import subprocess
import sys

BASE_ORDER = ['A', 'T', 'C', 'G', 'N', 'D']

CHROM_TO_CHR = {
    '2L': 'chr2L', '2R': 'chr2R',
    '3L': 'chr3L', '3R': 'chr3R',
    'X':  'chrX'
}
CHR_TO_CHROM = {v: k for k, v in CHROM_TO_CHR.items()}


def parse_markers(marker_file):
    """
    Read inversion_markers_v6.txt.
    Returns dict: {(chr_style_chrom, pos): inv_allele}
    e.g. {('chr2L', 1234567): 'A', ...}
    """
    markers = {}
    with open(marker_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            inv, chrom, pos, allele = parts[0], parts[1], int(parts[2]), parts[3].upper()
            chrom_chr = CHROM_TO_CHR.get(chrom, chrom)
            markers[(chrom_chr, pos)] = allele
    print(f"  {len(markers)} marker positions loaded from {marker_file}", file=sys.stderr)
    return markers


def get_ref_bases(fasta, positions):
    """
    Retrieve reference bases at all marker positions in one samtools faidx call.
    positions: list of (chr_style_chrom, pos) tuples
    Returns: dict {(chr_style_chrom, pos): base}
    """
    regions = [f"{chrom}:{pos}-{pos}" for chrom, pos in positions]
    result = subprocess.run(
        ['samtools', 'faidx', fasta] + regions,
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"ERROR: samtools faidx failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    ref_bases = {}
    current_key = None
    for line in result.stdout.split('\n'):
        if line.startswith('>'):
            # Header format: >chr2L:1234567-1234567
            header = line[1:]
            chrom_part, range_part = header.rsplit(':', 1)
            pos = int(range_part.split('-')[0])
            current_key = (chrom_part, pos)
        elif line.strip() and current_key:
            ref_bases[current_key] = line.strip().upper()[0]
            current_key = None

    print(f"  {len(ref_bases)} reference bases retrieved from {fasta}", file=sys.stderr)
    return ref_bases


def counts_to_sync(ref_count, alt_count, ref_base, inv_allele):
    """
    Build A:T:C:G:N:del string.
    - REF reads go into the ref_base slot.
    - ALT reads go into the inv_allele slot (if inv_allele != ref_base).
    - If inv_allele == ref_base (inversion allele is the reference), ALT reads go to N.
    """
    counts = {b: 0 for b in BASE_ORDER}
    if ref_base in counts:
        counts[ref_base] += ref_count
    if inv_allele != ref_base and inv_allele in counts:
        counts[inv_allele] += alt_count
    else:
        counts['N'] += alt_count
    return ':'.join(str(counts[b]) for b in BASE_ORDER)


def main():
    parser = argparse.ArgumentParser(description='Convert RefAlt to sync format for inversion typing.')
    parser.add_argument('--markers',    required=True, help='inversion_markers_v6.txt from DrosEU pipeline')
    parser.add_argument('--refalt-dir', required=True, help='Directory containing RefAlt.chrXX.txt files')
    parser.add_argument('--ref',        required=True, help='dm6 reference FASTA (must be samtools faidx indexed)')
    parser.add_argument('--output',     required=True, help='Output sync file path')
    parser.add_argument('--pops',       required=True, help='Output populations.txt file path (for InvFreq.py)')
    args = parser.parse_args()

    chroms_chr = ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX']

    # ---- Step 1: Load marker positions ----
    print("Loading marker positions...", file=sys.stderr)
    markers = parse_markers(args.markers)

    # ---- Step 2: Get reference bases at all marker positions ----
    print("Retrieving reference bases...", file=sys.stderr)
    ref_bases = get_ref_bases(args.ref, list(markers.keys()))

    # ---- Step 3: Stream RefAlt files, collect marker positions ----
    print("Reading RefAlt files...", file=sys.stderr)
    sample_names = None
    col_idx = {}
    refalt_data = {}  # {(chr_style_chrom, pos): {sample: (ref_c, alt_c)}}

    for chrom in chroms_chr:
        fname = os.path.join(args.refalt_dir, f"RefAlt.{chrom}.txt")
        if not os.path.exists(fname):
            print(f"  WARNING: {fname} not found, skipping", file=sys.stderr)
            continue

        print(f"  {fname}...", file=sys.stderr)
        with open(fname) as f:
            header = f.readline().strip().split()
            col_idx = {name: i for i, name in enumerate(header)}

            # Parse sample names from header columns on first file
            if sample_names is None:
                sample_names = [col[4:] for col in header if col.startswith('REF_')]
                print(f"  Samples detected: {sample_names}", file=sys.stderr)

            chrom_col = col_idx.get('CHROM', 0)
            pos_col   = col_idx.get('POS',   1)

            hits = 0
            for line in f:
                parts = line.strip().split()
                if not parts:
                    continue
                pos = int(parts[pos_col])
                row_chrom = parts[chrom_col]
                key = (row_chrom, pos)
                if key not in markers:
                    continue
                sample_counts = {}
                for samp in sample_names:
                    rc = int(parts[col_idx[f'REF_{samp}']])
                    ac = int(parts[col_idx[f'ALT_{samp}']])
                    sample_counts[samp] = (rc, ac)
                refalt_data[key] = sample_counts
                hits += 1
            print(f"    {hits} marker positions found on {chrom}", file=sys.stderr)

    n_found = len(refalt_data)
    n_total = len(markers)
    print(f"\nMarker positions found in RefAlt data: {n_found} / {n_total}", file=sys.stderr)
    if n_found < n_total * 0.5:
        print("WARNING: fewer than 50% of markers found. Check chromosome naming and coverage.",
              file=sys.stderr)

    # ---- Step 4: Write sync file ----
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    chrom_order = {c: i for i, c in enumerate(chroms_chr)}
    sorted_positions = sorted(
        refalt_data.keys(),
        key=lambda x: (chrom_order.get(x[0], 99), x[1])
    )

    print(f"\nWriting sync file: {args.output}", file=sys.stderr)
    with open(args.output, 'w') as out:
        for (chrom, pos) in sorted_positions:
            inv_allele = markers[(chrom, pos)]
            ref_base   = ref_bases.get((chrom, pos), 'N')
            chrom_out  = CHR_TO_CHROM.get(chrom, chrom)  # strip chr prefix for DrosEU

            sync_cols = [chrom_out, str(pos), ref_base]
            for samp in sample_names:
                rc, ac = refalt_data[(chrom, pos)].get(samp, (0, 0))
                sync_cols.append(counts_to_sync(rc, ac, ref_base, inv_allele))
            out.write('\t'.join(sync_cols) + '\n')

    print(f"Sync file written: {n_found} positions, {len(sample_names)} samples", file=sys.stderr)

    # ---- Step 5: Write populations.txt ----
    # InvFreq.py expects: population_id <tab> column_index (1-based, starting at 4)
    with open(args.pops, 'w') as pf:
        for i, samp in enumerate(sample_names):
            col = i + 4  # columns 1=chr, 2=pos, 3=ref, 4+=samples
            pf.write(f"{samp}\t{col}\n")
    print(f"Populations file written: {args.pops}", file=sys.stderr)
    print(f"Sample order (sync columns 4+): {sample_names}", file=sys.stderr)


if __name__ == '__main__':
    main()
