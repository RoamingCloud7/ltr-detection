
import sys
import os
import subprocess
import tempfile
import shutil
import argparse

import pandas as pd
import networkx as nx
import sourmash
from Bio import SeqIO

# -----------------------------------------------------------------------------
# I/O helpers
# -----------------------------------------------------------------------------
def load_candidates(tsv_path):
    """Load candidate TSV into a DataFrame."""
    return pd.read_csv(tsv_path, sep='\t')

def load_reads(fasta_path):
    """
    Load reads into a dict mapping only the first token of the header
    (everything before any whitespace) → SeqRecord.
    """
    reads = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        header_id = rec.id.split()[0]
        if header_id in reads:
            raise ValueError(f"Duplicate read ID in FASTA: {header_id}")
        reads[header_id] = rec
    return reads

# -----------------------------------------------------------------------------
# Anchor extraction
# -----------------------------------------------------------------------------
def extract_anchor_seqs(df, reads_dict):
    tmpdir = tempfile.mkdtemp(prefix="anchors_")
    fasta_paths = []
    for idx, row in df.iterrows():
        rid = row['seq_id']
        rec = reads_dict[rid]
        start, end = int(row['ltr1_end']), int(row['ltr2_start'])
        anchor = rec.seq[start:end]
        fa = os.path.join(tmpdir, f"anchor_{idx}.fa")
        with open(fa, 'w') as out:
            out.write(f">{idx}\n{anchor}\n")
        fasta_paths.append(fa)

    df2 = df.copy()
    df2['anchor_fasta'] = fasta_paths
    return df2, tmpdir

# -----------------------------------------------------------------------------
# MinHash sketching & bucketing
# -----------------------------------------------------------------------------
def compute_signatures(df, ksize=21, scaled=1000):
    """
    Compute a sourmash MinHash sketch for each anchor FASTA in df.
    Returns a list of MinHash objects aligned with df.index.
    """
    sigs = []
    for fpath in df['anchor_fasta']:
        mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
        for rec in SeqIO.parse(fpath, "fasta"):
            mh.add_sequence(str(rec.seq), force=True)
        sigs.append(mh)
    return sigs

def bucket_signatures(sigs, threshold=0.15):
    """
    Single-linkage bucketing of MinHash sketches: any pair with Jaccard ≤ thresh
    goes in the same bucket. Returns list of sets of indices.
    """
    n = len(sigs)
    buckets = []
    assigned = set()
    for i in range(n):
        if i in assigned:
            continue
        bucket = {i}
        for j in range(i+1, n):
            if j in assigned:
                continue
            if sigs[i].jaccard(sigs[j]) <= threshold:
                bucket.add(j)
        assigned.update(bucket)
        buckets.append(bucket)
    return buckets

# -----------------------------------------------------------------------------
# All-vs-all overlap & graph construction
# -----------------------------------------------------------------------------
def build_overlap_graph(bucket, df, preset="map-hifi"):
    """
    Concatenate anchor FASTAs in `bucket`, run minimap2 all-vs-all,
    and build an undirected graph where edges connect pairs with
    ≥95% identity over ≥80% of the shorter sequence.
    Nodes are integer indices from df.
    """
    G = nx.Graph()
    for idx in bucket:
        G.add_node(idx)

    # write a combined FASTA:
    tmpfa = tempfile.NamedTemporaryFile(delete=False, suffix=".fa").name
    with open(tmpfa, 'w') as out:
        for idx in bucket:
            with open(df.at[idx, 'anchor_fasta']) as inp:
                out.write(inp.read())

    # run minimap2
    cmd = ["minimap2", "-c", "--cs", "-x", preset, tmpfa, tmpfa]
    paf = subprocess.check_output(cmd).decode()

    # parse PAF
    for line in paf.splitlines():
        cols = line.split('\t')
        qname, qlen = cols[0], int(cols[1])
        tname, tlen = cols[5], int(cols[6])
        if qname == tname:
            continue
        aln_len = int(cols[3])
        pid = float(cols[9])
        if pid >= 95 and aln_len >= 0.8 * min(qlen, tlen):
            G.add_edge(int(qname), int(tname))

    os.unlink(tmpfa)
    return G

# -----------------------------------------------------------------------------
# Clustering & representative picking
# -----------------------------------------------------------------------------
def cluster_and_pick(df, sigs, min_support, preset):
    """
    - Bucket by MinHash
    - For each bucket, build overlap graph
    - Prune edges that are not part of any triangle
    - Extract connected components of pruned graph
    - Filter out components smaller than min_support
    - Pick highest-score member per component as representative
    - Return DataFrame of reps with support column
    """
    rep_idxs = []
    support = []

    buckets = bucket_signatures(sigs)
    for bucket in buckets:
        G = build_overlap_graph(bucket, df, preset=preset)

        # Prune edges that lack a common neighbor (break transitive chains)
        to_remove = []
        for u, v in G.edges():
            # require at least one common neighbor
            if len(set(G[u]).intersection(G[v])) < 1:
                to_remove.append((u, v))
        G.remove_edges_from(to_remove)

        # Extract clusters and pick reps
        for comp in nx.connected_components(G):
            if len(comp) < min_support:
                continue
            best = max(comp, key=lambda i: df.at[i, 'score'])
            rep_idxs.append(best)
            support.append(len(comp))

    out = df.loc[rep_idxs].copy()
    out['support'] = support
    return out

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Cluster per-read LTR-RT candidates and pick non-redundant reps."
    )
    p.add_argument("candidates_tsv",
                   help="Input TSV from candidate_in_reads")
    p.add_argument("reads_fasta",
                   help="Raw reads FASTA (CCS/HiFi or ONT)")
    p.add_argument("output_tsv",
                   help="Output TSV of representatives + support")
    p.add_argument("--min-support", type=int, default=3,
                   help="Discard clusters with support < this value")
    p.add_argument("--minhash-threshold", type=float, default=0.15,
                   help="Jaccard distance cutoff for MinHash bucketing")
    p.add_argument("--preset", default="map-hifi",
                   help="minimap2 preset for overlap (e.g. map-hifi, map-ont)")
    args = p.parse_args()

    # Load inputs
    df = load_candidates(args.candidates_tsv)
    reads = load_reads(args.reads_fasta)

    # Extract anchors
    df, tmpdir = extract_anchor_seqs(df, reads)

    # Sketch with MinHash
    sigs = compute_signatures(df)

    # Cluster + pick representatives
    reps = cluster_and_pick(df, sigs,
                            min_support=args.min_support,
                            preset=args.preset)

    # Filter and write output
    reps = reps[reps['support'] >= args.min_support]
    reps = reps.drop(columns=['anchor_fasta'])
    reps.to_csv(args.output_tsv, sep='\t', index=False)

    # Cleanup
    shutil.rmtree(tmpdir)
    print(f"Wrote {len(reps)} reps (support ≥ {args.min_support}) to {args.output_tsv}")

if __name__ == "__main__":
    main()
