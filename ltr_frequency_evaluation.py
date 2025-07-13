import pandas as pd
import os
import sys
from collections import Counter, defaultdict

def read_ref_bed_with_id(file_path):
    entries = []
    with open(file_path) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                chrom, start, end, id_ = line.split()[:4]
                entries.append((chrom, int(start), int(end), id_))
    return entries

def read_pred_bed(file_path):
    intervals = []
    with open(file_path) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                chrom, start, end = line.split()[:3]
                intervals.append((chrom, int(start), int(end)))
    return intervals

def merge_intervals(intervals):
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

def intersection_length(a_intervals, b_intervals):
    a_by_chr = defaultdict(list)
    b_by_chr = defaultdict(list)
    for chrom, s, e in a_intervals:
        a_by_chr[chrom].append((s, e))
    for chrom, s, e in b_intervals:
        b_by_chr[chrom].append((s, e))
    total = 0
    for chrom in a_by_chr:
        if chrom not in b_by_chr:
            continue
        a_list = sorted(a_by_chr[chrom])
        b_list = sorted(b_by_chr[chrom])
        i = j = 0
        while i < len(a_list) and j < len(b_list):
            a_s, a_e = a_list[i]
            b_s, b_e = b_list[j]
            lo, hi = max(a_s, b_s), min(a_e, b_e)
            if lo < hi:
                total += hi - lo
            if a_e < b_e:
                i += 1
            else:
                j += 1
    return total

if __name__ == "__main__":
    # Paths to your files
    ref_file = '/datadrive/maize/masking/Zm-Mo17-REFERENCE_msRepDB/maize_msRepDB.bed'
    pred_file = '/datadrive/maize/masking/Zm-Mo17-REFERENCE_repeatmodeler/maize_repeatmodeler.bed'

    # Read data
    ref_entries = read_ref_bed_with_id(ref_file)
    pred_intervals = read_pred_bed(pred_file)
    pred_merged = merge_intervals(pred_intervals)

    # Count occurrence frequency of each reference ID
    id_counts = Counter([id_ for _, _, _, id_ in ref_entries])
    counts = list(id_counts.values())

    # Statistics of frequency distribution
    stats = pd.Series(counts).describe()
    print("=== LTR-RT Copy Number Distribution ===")
    print(f"Total distinct IDs: {int(stats['count'])}")
    print(stats[['min', '25%', '50%', '75%', 'max', 'mean', 'std']].to_string())
    print("\nAssigning bins based on quartiles...")

    # Sort IDs by count, assign to high/med/low based on quartile
    sorted_ids = [id_ for id_, _ in id_counts.most_common()]
    total_ids = len(sorted_ids)
    H = total_ids // 4
    high_ids = set(sorted_ids[:H])
    low_ids = set(sorted_ids[-H:])
    med_ids = set(sorted_ids[H:-H])

    # Compute recall for each bin
    results = []
    for label, id_set in [('High-copy', high_ids), ('Medium-copy', med_ids), ('Low-copy', low_ids)]:
        # Filter reference entries by bin
        bin_entries = [(chrom, s, e) for chrom, s, e, id_ in ref_entries if id_ in id_set]
        if not bin_entries:
            continue
        bin_merged = merge_intervals(bin_entries)
        total_ref_bases = sum(e - s for _, s, e in bin_merged)
        tp_bases = intersection_length(pred_merged, bin_merged)
        recall = tp_bases / total_ref_bases if total_ref_bases else 0
        results.append({
            'Copy-frequency': label,
            'Num_IDs': len(id_set),
            'Total_ref_bases': total_ref_bases,
            'TP_bases': tp_bases,
            'Recall': recall
        })

    # Display recall results
    df = pd.DataFrame(results, columns=['Copy-frequency', 'Num_IDs', 'Total_ref_bases', 'TP_bases', 'Recall'])
    print("\n=== Stratified Recall by Copy-Frequency ===")
    print(df.to_string(index=False))
