import csv
import sys
from collections import defaultdict, namedtuple

Candidate = namedtuple("Candidate", [
    "seq_id", "ltr1_start", "ltr1_end",
    "ltr2_start", "ltr2_end",
    "identity", "score", "gap",
    "tsd5", "tsd3"
])

def overlaps(a_start, a_end, b_start, b_end):
    return not (a_end <= b_start or b_end <= a_start)

def filter_group(cands):
    """
    Given a list of Candidate for one seq_id, return a subset
    such that none of the kept candidates overlap each other
    on either the LTR1 interval or the LTR2 interval.
    We sort by descending score (then identity) and greedily select.
    """
    # sort by score desc, then identity desc
    cands_sorted = sorted(cands, key=lambda c: (-c.score, -c.identity))
    kept = []
    for c in cands_sorted:
        conflict = False
        for k in kept:
            if overlaps(c.ltr1_start, c.ltr1_end, k.ltr1_start, k.ltr1_end) \
               or overlaps(c.ltr2_start, c.ltr2_end, k.ltr2_start, k.ltr2_end):
                conflict = True
                break
        if not conflict:
            kept.append(c)
    # restore original order by ltr1_start if desired
    kept.sort(key=lambda c: c.ltr1_start)
    return kept

def main(input_tsv, output_tsv):
    groups = defaultdict(list)

    # Read all candidates
    with open(input_tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            c = Candidate(
                seq_id      = row["seq_id"],
                ltr1_start  = int(row["ltr1_start"]),
                ltr1_end    = int(row["ltr1_end"]),
                ltr2_start  = int(row["ltr2_start"]),
                ltr2_end    = int(row["ltr2_end"]),
                identity    = float(row["identity"]),
                score       = float(row["score"]),
                gap         = int(row["gap"]),
                tsd5        = row.get("tsd5",""),
                tsd3        = row.get("tsd3","")
            )
            groups[c.seq_id].append(c)

    # Filter each group
    with open(output_tsv, "w", newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        # write header
        writer.writerow([
            "seq_id","ltr1_start","ltr1_end",
            "ltr2_start","ltr2_end",
            "identity","score","gap",
            "tsd5","tsd3"
        ])
        for seq_id, cands in groups.items():
            kept = filter_group(cands)
            for c in kept:
                writer.writerow([
                    c.seq_id, c.ltr1_start, c.ltr1_end,
                    c.ltr2_start, c.ltr2_end,
                    f"{c.identity:.3f}", f"{c.score:.1f}", c.gap,
                    c.tsd5, c.tsd3
                ])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: filter_overlaps.py <candidates.tsv> <filtered.tsv>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])