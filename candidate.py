from Bio import SeqIO, pairwise2
from collections import namedtuple

# ------------------------------------------------------------------------
# Data structure to hold each candidate
# ------------------------------------------------------------------------
Candidate = namedtuple(
    "Candidate",
    [
        "seq_id",
        "ltr1_start",
        "ltr1_end",
        "ltr2_start",
        "ltr2_end",
        "identity",
        "score",
        "gap"
    ]
)

# ------------------------------------------------------------------------
# 1. k‐mer seeding
# ------------------------------------------------------------------------
def index_kmers(seq: str, k: int) -> dict:
    """
    Build a dict mapping each k‐mer in `seq` to a list of its positions.
    """
    kmer_dict = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        kmer_dict.setdefault(kmer, []).append(i)
    return kmer_dict

def find_seed_pairs(kmer_dict: dict, min_occ: int = 2) -> list:
    """
    From the kmer→positions map, emit all pairs of positions
    for k‐mers that occur >= min_occ times.
    """
    seeds = []
    for positions in kmer_dict.values():
        if len(positions) < min_occ:
            continue
        # every unique pair (i,j)
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                seeds.append((positions[i], positions[j]))
    return seeds

# ------------------------------------------------------------------------
# 2. cluster adjacent seed hits
# ------------------------------------------------------------------------
def cluster_seeds(
    seed_pairs: list,
    seed_len: int,
    max_gap: int
) -> list:
    """
    Merge any seed‐pairs whose windows (pos, pos+seed_len)
    overlap or lie within max_gap of each other.
    Returns a list of clusters: [ ((a1,a2),(b1,b2)), ... ]
    """
    # build initial windows
    intervals = [ (p1, p1+seed_len, p2, p2+seed_len)
                  for p1, p2 in seed_pairs ]
    intervals.sort(key=lambda x: (x[0], x[2]))
    used = [False]*len(intervals)
    clusters = []

    for i, (a1,a2,b1,b2) in enumerate(intervals):
        if used[i]:
            continue
        cur_a1, cur_a2 = a1, a2
        cur_b1, cur_b2 = b1, b2
        used[i] = True

        # try to absorb any nearby/overlapping intervals
        for j in range(i+1, len(intervals)):
            if used[j]:
                continue
            a1j,a2j,b1j,b2j = intervals[j]
            if (a1j <= cur_a2 + max_gap) and (b1j <= cur_b2 + max_gap):
                # merge
                cur_a1 = min(cur_a1, a1j)
                cur_a2 = max(cur_a2, a2j)
                cur_b1 = min(cur_b1, b1j)
                cur_b2 = max(cur_b2, b2j)
                used[j] = True

        clusters.append(((cur_a1,cur_a2),(cur_b1,cur_b2)))

    return clusters

# ------------------------------------------------------------------------
# 3. refinement via local alignment
# ------------------------------------------------------------------------
def refine_clusters(
    seq: str,
    clusters: list,
    buffer: int,
    match: float,
    mismatch: float,
    open_gap: float,
    extend_gap: float,
    min_identity: float
) -> list:
    """
    For each cluster window, do a Smith–Waterman alignment
    on the two buffered sub‐sequences. Keep only those
    with identity >= min_identity.
    Returns [ ((a1,a2),(b1,b2), identity, score), ... ]
    """
    refined = []
    for (a1,a2),(b1,b2) in clusters:
        w1 = seq[max(0,a1-buffer): min(len(seq),a2+buffer)]
        w2 = seq[max(0,b1-buffer): min(len(seq),b2+buffer)]
        aligns = pairwise2.align.localms(
            w1, w2,
            match, mismatch,
            open_gap, extend_gap,
            one_alignment_only=True
        )
        if not aligns:
            continue

        aln = aligns[0]
        seqA, seqB, score, _, _ = aln

        # compute percent identity over aligned columns
        matches = sum(
            1
            for x,y in zip(seqA, seqB)
            if x == y and x != "-" and y != "-"
        )
        length = sum(
            1
            for x,y in zip(seqA, seqB)
            if x != "-" and y != "-"
        )
        identity = (matches / length) if length > 0 else 0.0

        if identity >= min_identity:
            refined.append(((a1,a2),(b1,b2), identity, score))

    return refined

# ------------------------------------------------------------------------
# 4. distance & size checks
# ------------------------------------------------------------------------
def filter_refined(
    seq_id: str,
    refined: list,
    min_gap: int,
    max_gap: int,
    min_ltr: int,
    max_ltr: int,
    max_len_diff: float
) -> list:
    """
    From refined clusters, enforce:
      - LTR size in [min_ltr, max_ltr]
      - size difference <= max_len_diff fraction
      - inter‐LTR gap in [min_gap, max_gap]
      - (implicit) same orientation
    Returns a list of Candidate tuples.
    """
    cands = []
    for (a1,a2),(b1,b2), identity, score in refined:
        l1 = a2 - a1
        l2 = b2 - b1
        if not (min_ltr <= l1 <= max_ltr and min_ltr <= l2 <= max_ltr):
            continue
        if abs(l1 - l2) / min(l1, l2) > max_len_diff:
            continue
        gap = b1 - a2
        if gap < min_gap or gap > max_gap:
            continue

        cands.append(
            Candidate(
                seq_id=seq_id,
                ltr1_start=a1,
                ltr1_end=a2,
                ltr2_start=b1,
                ltr2_end=b2,
                identity=identity,
                score=score,
                gap=gap
            )
        )
    return cands

# ------------------------------------------------------------------------
# Driver function
# ------------------------------------------------------------------------
def identify_ltr_pairs(
    fasta_path: str,
    k: int = 5,
    min_identity: float = 0.6,
    min_gap: int = 100,
    max_gap: int = 15000,
    min_ltr: int = 100,
    max_ltr: int = 2000,
    max_len_diff: float = 0.10,
    buffer: int = 100,
    max_cluster_gap: int = 50,
    # alignment scoring
    match: float = 2,
    mismatch: float = -1,
    open_gap: float = -0.5,
    extend_gap: float = -0.1
) -> list:
    """
    Run the full LTR‐pair candidate identification on every record
    in the input FASTA. Returns a flat list of Candidate objects.
    """
    all_candidates = []

    count = 0

    for rec in SeqIO.parse(fasta_path, "fasta"):
        count = count + 1
        print("analyse the rec " + str(count))
        seq = str(rec.seq)
        # 1) seed
        kmers   = index_kmers(seq, k)
        seeds   = find_seed_pairs(kmers)
        # 2) cluster
        clusters = cluster_seeds(seeds, k, max_cluster_gap)
        # 3) refine
        refined  = refine_clusters(
            seq, clusters, buffer,
            match, mismatch, open_gap, extend_gap,
            min_identity
        )
        # 4) filter
        filtered = filter_refined(
            rec.id, refined,
            min_gap, max_gap,
            min_ltr, max_ltr,
            max_len_diff
        )
        all_candidates.extend(filtered)

    return all_candidates

# ------------------------------------------------------------------------
# CLI interface
# ------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Identify LTR‐pair candidates in repetitive sequences"
    )
    parser.add_argument("-i", "--input",  required=True,
                        help="FASTA of repeats")
    parser.add_argument("-o", "--output", required=True,
                        help="TSV output of candidates")
    parser.add_argument("--k",            type=int,   default=5,
                        help="k‐mer size")
    parser.add_argument("--min_identity", type=float, default=0.6,
                        help="min % identity for LTR pair")
    parser.add_argument("--min_gap",      type=int,   default=100,
                        help="min bp between LTRs")
    parser.add_argument("--max_gap",      type=int,   default=15000,
                        help="max bp between LTRs")
    parser.add_argument("--min_ltr",      type=int,   default= 100,
                        help="min LTR length")
    parser.add_argument("--max_ltr",      type=int,   default=2000,
                        help="max LTR length")
    parser.add_argument("--buffer",       type=int,   default=100,
                        help="flanking buffer for alignment")
    args = parser.parse_args()

    candidates = identify_ltr_pairs(
        fasta_path    = args.input,
        k              = args.k,
        min_identity   = args.min_identity,
        min_gap        = args.min_gap,
        max_gap        = args.max_gap,
        min_ltr        = args.min_ltr,
        max_ltr        = args.max_ltr,
        buffer         = args.buffer
    )

    # write TSV
    with open(args.output, "w") as out:
        out.write(
            "seq_id\tltr1_start\tltr1_end\tltr2_start\tltr2_end\t"
            "identity\tscore\tgap\n"
        )
        for c in candidates:
            out.write(
                f"{c.seq_id}\t{c.ltr1_start}\t{c.ltr1_end}\t"
                f"{c.ltr2_start}\t{c.ltr2_end}\t"
                f"{c.identity:.3f}\t{c.score:.1f}\t{c.gap}\n"
            )
    print(f"Wrote {len(candidates)} LTR‐pair candidates to {args.output}")


