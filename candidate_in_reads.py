import sys
import csv
import argparse
import itertools
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from Bio.Align import PairwiseAligner

Candidate = namedtuple(
    "Candidate",
    [
        "read_id", "ltr1_start", "ltr1_end",
        "ltr2_start", "ltr2_end",
        "identity", "score", "gap"
    ]
)

# ------------------------------------------------------------------------
# 1. k‐mer seeding
# ------------------------------------------------------------------------
def index_kmers(seq: str, k: int) -> dict:
    d = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        d.setdefault(kmer, []).append(i)
    return d

def find_seed_pairs(kmer_dict: dict, min_occ: int = 2) -> list:
    seeds = []
    for pos_list in kmer_dict.values():
        if len(pos_list) < min_occ:
            continue
        for i in range(len(pos_list)):
            for j in range(i+1, len(pos_list)):
                seeds.append((pos_list[i], pos_list[j]))
    return seeds

# ------------------------------------------------------------------------
# 2. cluster adjacent seeds
# ------------------------------------------------------------------------
def cluster_seeds(seed_pairs: list, seed_len: int, max_gap: int) -> list:
    intervals = [(a, a+seed_len, b, b+seed_len) for a,b in seed_pairs]
    intervals.sort(key=lambda x: (x[0], x[2]))
    used = [False]*len(intervals)
    clusters = []
    for i,(a1,a2,b1,b2) in enumerate(intervals):
        if used[i]: continue
        cur_a1,cur_a2,cur_b1,cur_b2 = a1,a2,b1,b2
        used[i] = True
        for j in range(i+1, len(intervals)):
            if used[j]: continue
            a1j,a2j,b1j,b2j = intervals[j]
            if a1j <= cur_a2 + max_gap and b1j <= cur_b2 + max_gap:
                cur_a1,cur_a2 = min(cur_a1,a1j), max(cur_a2,a2j)
                cur_b1,cur_b2 = min(cur_b1,b1j), max(cur_b2,b2j)
                used[j] = True
        clusters.append(((cur_a1,cur_a2),(cur_b1,cur_b2)))
    return clusters

# ------------------------------------------------------------------------
# 4. distance & size filtering
# ------------------------------------------------------------------------
def filter_refined(refined: list,
                   min_gap:int, max_gap:int,
                   min_ltr:int, max_ltr:int,
                   max_len_diff:float) -> list:
    out = []
    for (a1,a2),(b1,b2),identity,score in refined:
        l1, l2 = a2-a1, b2-b1
        if not (min_ltr<=l1<=max_ltr and min_ltr<=l2<=max_ltr):
            continue
        if abs(l1-l2)/min(l1,l2) > max_len_diff:
            continue
        gap = b1 - a2
        if gap<min_gap or gap>max_gap:
            continue
        out.append((a1,a2,b1,b2,identity,score,gap))
    return out

# ------------------------------------------------------------------------
# helper: find TSD 2–20 bp immediately flanking [start,end)
# ------------------------------------------------------------------------
def find_tsd(seq: str, start: int, end: int,
             min_len: int = 2, max_len: int = 20) -> tuple:
    for L in range(max_len, min_len-1, -1):
        if start-L<0 or end+L>len(seq): continue
        left  = seq[start-L:start]
        right = seq[end:end+L]
        if left == right:
            return left,right
    return "",""

# ------------------------------------------------------------------------
# per-read worker
# ------------------------------------------------------------------------
def process_record(rec, params):
    seq = str(rec.seq)
    L = len(seq)
    if L < params["min_read_len"] or L > params["max_read_len"]:
        return []
    # 1) seed
    seeds = find_seed_pairs(index_kmers(seq, params["k"]))
    # 2) cluster
    clusters = cluster_seeds(seeds, params["k"], params["max_cluster_gap"])
    # 3) refine via local alignment
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score   = params["match"]
    aligner.mismatch_score= params["mismatch"]
    aligner.open_gap_score= params["open_gap"]
    aligner.extend_gap_score = params["extend_gap"]

    refined = []
    for (a1,a2),(b1,b2) in clusters:
        w1 = seq[max(0,a1-params["buffer"]): a2+params["buffer"]]
        w2 = seq[max(0,b1-params["buffer"]): b2+params["buffer"]]
        it = aligner.align(w1,w2)
        try:
            aln = next(it)
        except StopIteration:
            continue
        # compute identity
        matches = 0
        aligned_len = 0
        segs1, segs2 = aln.aligned
        for (s1,e1),(s2,e2) in zip(segs1,segs2):
            length = e1 - s1
            aligned_len += length
            for o in range(length):
                if w1[s1+o] == w2[s2+o]:
                    matches += 1
        identity = matches/aligned_len if aligned_len else 0.0
        if identity >= params["min_identity"]:
            refined.append(((a1,a2),(b1,b2),identity,aln.score))

    # 4) filter
    filtered = filter_refined(refined,
                              params["min_gap"], params["max_gap"],
                              params["min_ltr"], params["max_ltr"],
                              params["max_len_diff"])
    # annotate TSD
    out = []
    for a1,a2,b1,b2,identity,score,gap in filtered:
        tsd5, tsd3 = find_tsd(seq, a1, b2)
        out.append((rec.id, a1, a2, b1, b2,
                    f"{identity:.3f}", f"{score:.1f}", gap,
                    tsd5, tsd3))
    return out

# ------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Parallel LTR‐pair detection in long reads"
    )
    p.add_argument("-i","--input", required=True,
                   help="FASTA/FASTQ of reads")
    p.add_argument("-o","--output",required=True,
                   help="TSV out")
    p.add_argument("--threads", type=int, default=32,
                   help="# parallel workers")
    p.add_argument("--min_read_len", type=int, default=1000)
    p.add_argument("--max_read_len", type=int, default=200000)
    p.add_argument("--k", type=int, default=11)
    p.add_argument("--min_identity", type=float, default=0.70)
    p.add_argument("--min_gap", type=int, default=100)
    p.add_argument("--max_gap", type=int, default=20000)
    p.add_argument("--min_ltr", type=int, default=100)
    p.add_argument("--max_ltr", type=int, default=2000)
    p.add_argument("--max_len_diff", type=float, default=0.20)
    p.add_argument("--buffer", type=int, default=100)
    p.add_argument("--max_cluster_gap", type=int, default=20)
    p.add_argument("--match", type=float, default=2.0)
    p.add_argument("--mismatch", type=float, default=-0.5)
    p.add_argument("--open_gap", type=float, default=-0.3)
    p.add_argument("--extend_gap", type=float, default=-0.05)
    args = p.parse_args()

    fmt = "fastq" if args.input.lower().endswith((".fq",".fastq")) else "fasta"

    params = {
        "min_read_len": args.min_read_len,
        "max_read_len": args.max_read_len,
        "k":             args.k,
        "min_identity":  args.min_identity,
        "min_gap":       args.min_gap,
        "max_gap":       args.max_gap,
        "min_ltr":       args.min_ltr,
        "max_ltr":       args.max_ltr,
        "max_len_diff":  args.max_len_diff,
        "buffer":        args.buffer,
        "max_cluster_gap":args.max_cluster_gap,
        "match":         args.match,
        "mismatch":      args.mismatch,
        "open_gap":      args.open_gap,
        "extend_gap":    args.extend_gap
    }

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "seq_id","ltr1_start","ltr1_end",
            "ltr2_start","ltr2_end",
            "identity","score","gap",
            "tsd5","tsd3"
        ])

        records = SeqIO.parse(args.input, fmt)
        with ProcessPoolExecutor(max_workers=args.threads) as exe:
            for result in exe.map(process_record,
                                  records,
                                  itertools.repeat(params)):
                for row in result:
                    writer.writerow(row)

    print("Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
