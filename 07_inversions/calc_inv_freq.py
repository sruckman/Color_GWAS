#!/usr/bin/env python3
"""
calc_inv_freq.py

Calculate inversion frequencies from a sync file filtered to marker positions.
Replaces InvFreq.py (Python 2) with a pure Python 3 implementation.

For each marker position, extracts the count of the inversion-associated allele
from the sync file and calculates frequency = inv_count / total_count.
Averages across all markers for each inversion x population combination.

Usage:
    python3 calc_inv_freq.py \
        --sync  process/inversions/all_populations.sync \
        --inv   process/inversions/inversion_markers_v6.txt \
        --pop   process/inversions/populations.txt \
        > process/inversions/inversion_frequencies.txt
"""

import argparse
import sys
from collections import defaultdict

BASE_ORDER = ['A', 'T', 'C', 'G', 'N', 'D']

def parse_sync_counts(field):
    """Parse A:T:C:G:N:del string into a dict."""
    vals = list(map(int, field.split(':')))
    return dict(zip(BASE_ORDER, vals))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sync', required=True, help='Sync file (marker positions only)')
    parser.add_argument('--inv',  required=True, help='inversion_markers_v6.txt')
    parser.add_argument('--pop',  required=True, help='populations.txt (name <tab> column, 1-based)')
    args = parser.parse_args()

    # Load populations: name -> 0-based index into sync count columns
    pops = []
    pop_col = {}
    with open(args.pop) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            name = parts[0]
            col  = int(parts[1]) - 1  # convert 1-based to 0-based
            pops.append(name)
            pop_col[name] = col

    # Load marker file: {(chrom, pos): (inversion_name, inv_allele)}
    markers = {}
    with open(args.inv) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            inv, chrom, pos, allele = parts[0], parts[1], int(parts[2]), parts[3].upper()
            markers[(chrom, pos)] = (inv, allele)

    # inversions list in order
    inv_names = []
    seen = set()
    with open(args.inv) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 1:
                continue
            inv = parts[0]
            if inv not in seen:
                inv_names.append(inv)
                seen.add(inv)

    # Accumulate per-inversion per-population allele frequencies
    # freq_sum[inv][pop] = list of per-marker frequencies
    freq_data = defaultdict(lambda: defaultdict(list))
    marker_counts = defaultdict(int)  # how many markers found per inversion

    with open(args.sync) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            chrom = fields[0]
            pos   = int(fields[1])
            # ref_base = fields[2]  # not needed here
            count_fields = fields[3:]

            key = (chrom, pos)
            if key not in markers:
                continue

            inv_name, inv_allele = markers[key]
            marker_counts[inv_name] += 1

            for pop in pops:
                col = pop_col[pop] - 3  # count_fields index: sync col - 3 (chr,pos,ref = cols 0-2)
                if col < 0 or col >= len(count_fields):
                    continue
                counts = parse_sync_counts(count_fields[col])
                inv_count = counts.get(inv_allele, 0)
                total     = sum(counts[b] for b in ['A', 'T', 'C', 'G'])
                if total == 0:
                    continue  # no coverage at this position
                freq = inv_count / total
                freq_data[inv_name][pop].append(freq)

    # Print results
    header = 'Inv\t' + '\t'.join(pops)
    print(header)

    for inv in inv_names:
        row = [inv]
        n_markers = marker_counts.get(inv, 0)
        for pop in pops:
            freqs = freq_data[inv][pop]
            if freqs:
                avg = sum(freqs) / len(freqs)
                row.append(f'{avg:.4f}')
            else:
                row.append('NA')
        print('\t'.join(row))
        print(f'  # {inv}: {n_markers} markers used', file=sys.stderr)

if __name__ == '__main__':
    main()
