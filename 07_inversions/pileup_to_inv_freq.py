#!/usr/bin/env python3
"""
pileup_to_inv_freq.py

Parse samtools mpileup output (multiple BAMs) at inversion marker positions
and calculate per-inversion frequencies per sample.

Mpileup columns per sample (after chr, pos, ref): depth, bases, quals
Base string codes:
    . ,         = reference allele (fwd/rev strand)
    A T C G     = specific alt base (fwd strand)
    a t c g     = specific alt base (rev strand)
    ^ $         = read start/end markers (skip)
    + - N       = indels (skip)
    *           = deletion placeholder (skip)

Usage:
    python3 pileup_to_inv_freq.py \\
        --pileup  process/inversions/marker_pileup.txt \\
        --markers process/inversions/inversion_markers_v6.txt \\
        --samples HOULE_L1F HOULE_L2F HOULE_L3F HOUSTON_L1F HOUSTON_L2F HOUSTON_L3F \\
        > process/inversions/inversion_frequencies_bam.txt
"""

import argparse
import re
import sys
from collections import defaultdict


def parse_base_string(base_str, ref_base):
    """
    Count A/T/C/G occurrences in a pileup base string.
    . and , are counted as the reference base.
    Indel sequences and read markers are skipped.
    Returns a dict {base: count}.
    """
    counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
    i = 0
    s = base_str.upper()
    ref = ref_base.upper()

    while i < len(s):
        c = s[i]
        if c in ('.', ','):
            if ref in counts:
                counts[ref] += 1
            i += 1
        elif c in counts:
            counts[c] += 1
            i += 1
        elif c == '^':
            i += 2  # skip ^ and the next mapping quality char
        elif c == '$':
            i += 1
        elif c in ('+', '-'):
            # skip indel: +Nbases or -Nbases
            i += 1
            num_str = ''
            while i < len(s) and s[i].isdigit():
                num_str += s[i]
                i += 1
            if num_str:
                i += int(num_str)  # skip the indel bases
        elif c == '*':
            i += 1  # deletion placeholder
        elif c == 'N':
            i += 1
        else:
            i += 1

    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pileup',  required=True, help='samtools mpileup output file')
    parser.add_argument('--markers', required=True, help='inversion_markers_v6.txt')
    parser.add_argument('--samples', required=True, nargs='+', help='Sample names in BAM order')
    args = parser.parse_args()

    samples = args.samples
    n_samples = len(samples)

    # Load marker file: {(chrom, pos): (inv_name, inv_allele)}
    # Marker file uses 2L/2R style; pileup uses chr2L style
    markers = {}
    inv_names = []
    seen_inv = set()
    with open(args.markers) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            inv, chrom, pos, allele = parts[0], parts[1], int(parts[2]), parts[3].upper()
            chrom_chr = 'chr' + chrom  # convert 2L -> chr2L
            markers[(chrom_chr, pos)] = (inv, allele)
            if inv not in seen_inv:
                inv_names.append(inv)
                seen_inv.add(inv)

    print(f"Loaded {len(markers)} marker positions for {len(inv_names)} inversions",
          file=sys.stderr)

    # Accumulate per-inversion per-sample allele frequencies and raw counts
    freq_data    = defaultdict(lambda: defaultdict(list))
    inv_counts   = defaultdict(lambda: defaultdict(int))  # total inv allele reads
    total_counts = defaultdict(lambda: defaultdict(int))  # total reads
    marker_used  = defaultdict(int)
    marker_miss  = defaultdict(int)

    with open(args.pileup) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')

            chrom  = fields[0]
            pos    = int(fields[1])
            ref    = fields[2].upper()

            key = (chrom, pos)
            if key not in markers:
                continue

            inv_name, inv_allele = markers[key]

            # Fields after chr/pos/ref: groups of 3 per sample (depth, bases, quals)
            sample_fields = fields[3:]
            if len(sample_fields) < n_samples * 3:
                marker_miss[inv_name] += 1
                continue

            any_coverage = False
            for i, samp in enumerate(samples):
                depth_str = sample_fields[i * 3]
                base_str  = sample_fields[i * 3 + 1]
                depth = int(depth_str)
                if depth == 0:
                    continue

                counts = parse_base_string(base_str, ref)
                total  = sum(counts.values())
                if total == 0:
                    continue

                inv_count = counts.get(inv_allele, 0)
                freq = inv_count / total
                freq_data[inv_name][samp].append(freq)
                inv_counts[inv_name][samp]   += inv_count
                total_counts[inv_name][samp] += total
                any_coverage = True

            if any_coverage:
                marker_used[inv_name] += 1
            else:
                marker_miss[inv_name] += 1

    # Print summary to stderr
    print("\nMarkers used per inversion:", file=sys.stderr)
    for inv in inv_names:
        print(f"  {inv}: {marker_used.get(inv, 0)} used, "
              f"{marker_miss.get(inv, 0)} missing/no coverage", file=sys.stderr)

    # Print frequency table
    header = 'Inv\t' + '\t'.join(samples)
    print(header)

    for inv in inv_names:
        row = [inv]
        for samp in samples:
            freqs = freq_data[inv][samp]
            if freqs:
                avg = sum(freqs) / len(freqs)
                row.append(f'{avg:.4f}')
            else:
                row.append('NA')
        print('\t'.join(row))

    # Write counts table (for Fisher's exact tests)
    counts_file = args.pileup.replace('.txt', '_counts.txt')
    with open(counts_file, 'w') as cf:
        cf.write('Inv\tSample\tinv_reads\ttotal_reads\n')
        for inv in inv_names:
            for samp in samples:
                ic = inv_counts[inv][samp]
                tc = total_counts[inv][samp]
                cf.write(f'{inv}\t{samp}\t{ic}\t{tc}\n')
    print(f"\nCounts: {counts_file}", file=sys.stderr)


if __name__ == '__main__':
    main()
