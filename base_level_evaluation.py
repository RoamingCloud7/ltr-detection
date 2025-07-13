import pandas as pd

import sys
from collections import defaultdict

def read_bed(file_path):
    """Read a BED file into a list of (chrom, start, end) tuples."""
    intervals = []
    with open(file_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            chrom, start, end = line.split()[:3]
            intervals.append((chrom, int(start), int(end)))
    return intervals

def merge_intervals(intervals):
    """
    Merge overlapping or adjacent intervals _per chromosome_.
    Returns a new list of merged (chrom, start, end).
    """
    by_chr = defaultdict(list)
    for chrom, s, e in intervals:
        by_chr[chrom].append((s, e))

    merged = []
    for chrom, ivs in by_chr.items():
        ivs.sort()
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s <= cur_e:        
                cur_e = max(cur_e, e)
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((chrom, cur_s, cur_e))
    return merged

def total_length(intervals):
    return sum(e - s for _, s, e in intervals)

def intersection_length(a_intervals, b_intervals):
    """
    Two-pointer scan per chromosome to get total overlapping bases.
    """
    a_by_chr = defaultdict(list)
    b_by_chr = defaultdict(list)
    for chrom, s, e in a_intervals:
        a_by_chr[chrom].append((s, e))
    for chrom, s, e in b_intervals:
        b_by_chr[chrom].append((s, e))

    total = 0
    for chrom, a_list in a_by_chr.items():
        if chrom not in b_by_chr:
            continue
        b_list = b_by_chr[chrom]
        a_list.sort(); b_list.sort()
        i = j = 0
        while i < len(a_list) and j < len(b_list):
            a_s, a_e = a_list[i]
            b_s, b_e = b_list[j]
            # overlap?
            lo = max(a_s, b_s)
            hi = min(a_e, b_e)
            if lo < hi:
                total += hi - lo
            # advance the one that ends first
            if a_e < b_e:
                i += 1
            else:
                j += 1
    return total

if __name__ == "__main__":
    # Paths to merged BED files
    file_a = '/datadrive/maize/masking/Zm-Mo17-REFERENCE_repeatmodeler/maize_repeatmodeler.bed'
    file_b = '/datadrive/maize/masking/Zm-Mo17-REFERENCE_msRepDB/maize_msRepDB.bed'

    # Read & merge
    a_raw = read_bed(file_a)
    b_raw = read_bed(file_b)
    a_ints = merge_intervals(a_raw)
    b_ints = merge_intervals(b_raw)

    # Base‐level counts
    sumA = total_length(a_ints)         # TP + FP
    sumB = total_length(b_ints)         # TP + FN
    TP   = intersection_length(a_ints, b_ints)

    # Clamp to avoid negatives
    FP = max(0, sumA - TP)
    FN = max(0, sumB - TP)

    # Metrics
    precision = TP / (TP + FP) if TP + FP > 0 else 0
    recall    = TP / (TP + FN) if TP + FN > 0 else 0
    f1        = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    # Report
    print(f"Total masked A (sumA): {sumA}")
    print(f"Total masked B (sumB): {sumB}")
    print(f"TP bases:              {TP}")
    print(f"FP bases:              {FP}")
    print(f"FN bases:              {FN}\n")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1 score:  {f1:.4f}")