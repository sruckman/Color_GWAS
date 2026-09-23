#!/usr/bin/env python3
"""
remap_markers.py

Sequence-based remapping of Kapun 2014 inversion marker SNPs from dm3 to dm6.

For each marker:
  1. Extracts a ±50bp window from the dm3 reference using samtools faidx
  2. Aligns that sequence to dm6 using BWA MEM
  3. Parses the SAM output to determine the dm6 coordinate of the center base
  4. Verifies the base at the dm6 coordinate matches the expected marker allele
  5. Compares to liftOver result (if available)

Outputs:
  - inversion_markers_v6_remapped.txt  : corrected marker file
  - remap_report.txt                   : per-marker comparison table

Usage:
    python3 remap_markers.py \\
        --markers  process/inversions/inversion_markers_v6.txt \\
        --dm3      ref/dm3.fa \\
        --dm6      ref/dm6.fa \\
        --out-dir  process/inversions
"""

import argparse
import os
import subprocess
import sys
import tempfile

FLANK = 50  # bp on each side of marker position


def get_seq(fasta, chrom, start, end):
    """
    Extract sequence from fasta using samtools faidx.
    start/end are 1-based inclusive.
    Returns uppercase sequence string or None on failure.
    """
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        ['samtools', 'faidx', fasta, region],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return None
    lines = result.stdout.strip().split('\n')
    if len(lines) < 2:
        return None
    return ''.join(lines[1:]).upper()


def bwa_mem_align(query_fa, ref_fa, bwa_idx_prefix=None):
    """
    Align sequences in query_fa against ref_fa using bwa mem.
    Returns list of SAM records (non-header lines).
    """
    target = bwa_idx_prefix if bwa_idx_prefix else ref_fa
    result = subprocess.run(
        ['bwa', 'mem', '-k', '19', '-T', '50', target, query_fa],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"BWA error: {result.stderr[:200]}", file=sys.stderr)
        return []
    return [l for l in result.stdout.split('\n')
            if l and not l.startswith('@')]


def parse_cigar_offset(cigar, query_pos):
    """
    Given a CIGAR string and a 0-based position in the query,
    return the corresponding offset in the reference.
    Returns None if the query position is in a deletion or clipped region.
    """
    import re
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    ref_offset = 0
    qry_offset = 0
    for length, op in ops:
        length = int(length)
        if op in ('M', '=', 'X'):
            if qry_offset + length > query_pos:
                ref_offset += (query_pos - qry_offset)
                return ref_offset
            qry_offset += length
            ref_offset += length
        elif op == 'I':
            if qry_offset + length > query_pos:
                return None  # insertion — no ref position
            qry_offset += length
        elif op in ('D', 'N'):
            ref_offset += length
        elif op in ('S', 'H'):
            if op == 'S':
                if qry_offset + length > query_pos:
                    return None  # soft-clipped
                qry_offset += length
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--markers',  required=True,
                        help='Marker file in dm3 coordinates (tab: inv chr pos allele)')
    parser.add_argument('--dm3',      required=True,
                        help='dm3 reference FASTA (samtools faidx indexed)')
    parser.add_argument('--dm6',      required=True,
                        help='dm6 reference FASTA (BWA indexed)')
    parser.add_argument('--liftover', default=None,
                        help='Optional: liftOver result file for comparison')
    parser.add_argument('--out-dir',  default='.',
                        help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # dm3 uses chr-style names (chr2L etc.)
    CHROM_MAP = {
        '2L': 'chr2L', '2R': 'chr2R',
        '3L': 'chr3L', '3R': 'chr3R',
        'X':  'chrX'
    }

    # Load markers (dm3 coordinates)
    markers = []
    with open(args.markers) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            inv, chrom, pos, allele = parts[0], parts[1], int(parts[2]), parts[3].upper()
            chrom_chr = CHROM_MAP.get(chrom, chrom)
            markers.append((inv, chrom, chrom_chr, pos, allele))

    print(f"Processing {len(markers)} markers...", file=sys.stderr)

    # Load liftOver results if provided
    liftover_coords = {}
    if args.liftover and os.path.exists(args.liftover):
        with open(args.liftover) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    name = parts[3]  # inv|chr|allele
                    lo_chrom = parts[0]
                    lo_pos = int(parts[2])  # BED end = 1-based pos
                    liftover_coords[name] = (lo_chrom, lo_pos)

    # Extract sequences from dm3 and write query FASTA
    print("Extracting dm3 sequences...", file=sys.stderr)
    query_seqs = {}  # {marker_id: (seq, center_in_seq)}

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa',
                                     delete=False) as tmpfa:
        query_fa = tmpfa.name
        for inv, chrom, chrom_chr, pos, allele in markers:
            marker_id = f"{inv}|{chrom}|{pos}|{allele}"
            start = max(1, pos - FLANK)
            end   = pos + FLANK

            seq = get_seq(args.dm3, chrom_chr, start, end)
            if seq is None:
                print(f"  WARNING: could not extract {chrom_chr}:{start}-{end}",
                      file=sys.stderr)
                continue

            center = pos - start  # 0-based index of marker base in seq
            query_seqs[marker_id] = (seq, center, pos)

            # Write FASTA entry — name encodes marker_id and center position
            safe_id = marker_id.replace('(', '_').replace(')', '_')
            tmpfa.write(f">{safe_id}|center={center}\n{seq}\n")

    # Align against dm6 with BWA
    print("Aligning to dm6 with BWA MEM...", file=sys.stderr)
    sam_records = bwa_mem_align(query_fa, args.dm6)
    os.unlink(query_fa)

    # Parse SAM to get dm6 coordinates
    dm6_coords = {}  # {marker_id: (chrom, pos, mapped_allele)}

    for rec in sam_records:
        fields = rec.split('\t')
        if len(fields) < 10:
            continue
        qname  = fields[0]
        flag   = int(fields[1])
        rname  = fields[2]
        pos_1  = int(fields[3])  # 1-based leftmost mapping position
        cigar  = fields[5]
        seq    = fields[9]

        if flag & 4:  # unmapped
            continue
        if rname == '*':
            continue

        # Parse marker_id and center from qname
        # Format: inv_chr_pos_allele|center=N
        try:
            parts = qname.rsplit('|center=', 1)
            marker_part = parts[0]
            center_in_query = int(parts[1])
        except (IndexError, ValueError):
            continue

        # Reverse complement affects center position
        if flag & 16:  # reverse strand
            seq_len = len(seq)
            center_in_query = seq_len - 1 - center_in_query

        # Get reference offset for center base
        ref_offset = parse_cigar_offset(cigar, center_in_query)
        if ref_offset is None:
            continue

        dm6_pos = pos_1 + ref_offset  # 1-based

        # Get base at dm6 position from reference
        dm6_base = get_seq(args.dm6, rname, dm6_pos, dm6_pos)

        # Recover original marker_id (undo safe replacement)
        # Try to match back to original markers
        orig_id = None
        for inv, chrom, chrom_chr, pos, allele in markers:
            safe = f"{inv}|{chrom}|{pos}|{allele}".replace('(', '_').replace(')', '_')
            if safe == marker_part:
                orig_id = f"{inv}|{chrom}|{pos}|{allele}"
                break

        if orig_id:
            dm6_coords[orig_id] = (rname, dm6_pos, dm6_base or '?')

    # Write report and corrected marker file
    report_path = os.path.join(args.out_dir, 'remap_report.txt')
    out_path    = os.path.join(args.out_dir, 'inversion_markers_v6_remapped.txt')
    pos_path    = os.path.join(args.out_dir, 'inversion_markers_v6_remapped_pos.txt')

    n_match = 0
    n_mismatch = 0
    n_unmapped = 0

    CHROM_STRIP = {v: k for k, v in CHROM_MAP.items()}

    with open(report_path, 'w') as rep, \
         open(out_path, 'w') as out, \
         open(pos_path, 'w') as pos_out:

        rep.write('\t'.join(['inv', 'chrom', 'dm3_pos', 'allele',
                             'dm6_chrom', 'dm6_pos', 'dm6_ref_base',
                             'allele_match', 'liftover_pos', 'liftover_agrees']) + '\n')

        for inv, chrom, chrom_chr, pos, allele in markers:
            marker_id = f"{inv}|{chrom}|{pos}|{allele}"

            lo_key = marker_id.replace('(', '_').replace(')', '_') if args.liftover else None

            if marker_id in dm6_coords:
                dm6_chrom, dm6_pos, dm6_base = dm6_coords[marker_id]
                # allele match: marker allele could be REF or ALT
                allele_match = (dm6_base == allele)

                lo_pos = liftover_coords.get(lo_key, (None, None))[1]
                lo_agrees = (lo_pos == dm6_pos) if lo_pos else 'N/A'

                if allele_match:
                    n_match += 1
                else:
                    n_mismatch += 1

                rep.write('\t'.join([inv, chrom, str(pos), allele,
                                     dm6_chrom, str(dm6_pos), dm6_base or '?',
                                     str(allele_match), str(lo_pos or 'N/A'),
                                     str(lo_agrees)]) + '\n')

                # Write corrected marker file
                chrom_short = CHROM_STRIP.get(dm6_chrom, dm6_chrom)
                out.write(f"{inv}\t{chrom_short}\t{dm6_pos}\t{allele}\n")
                pos_out.write(f"{chrom_short}\t{dm6_pos}\n")

            else:
                n_unmapped += 1
                rep.write('\t'.join([inv, chrom, str(pos), allele,
                                     'UNMAPPED', 'NA', 'NA',
                                     'NA', 'NA', 'NA']) + '\n')

    print(f"\nResults:", file=sys.stderr)
    print(f"  Mapped:        {len(dm6_coords)}", file=sys.stderr)
    print(f"  Allele match:  {n_match}", file=sys.stderr)
    print(f"  Allele mismatch (opposite strand or wrong pos): {n_mismatch}", file=sys.stderr)
    print(f"  Unmapped:      {n_unmapped}", file=sys.stderr)
    print(f"\nReport:  {report_path}", file=sys.stderr)
    print(f"Markers: {out_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
