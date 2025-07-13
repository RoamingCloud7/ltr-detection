#!/usr/bin/env python3
"""
truncated_verification.py

Phase II: Scan leftover repetitive regions (without high‐confidence LTR pairs)
for partial (truncated) LTR-RT signatures, using the output of your LTR-pair
detector to exclude full-length calls.

Produces:
 - TSV of truncated candidates with scores and coordinates.
 - FASTA of truncated regions with headers:
     "<seq_id>":<start>_<end>
"""

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path
from collections import namedtuple

from Bio import SeqIO
from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord

# -----------------------------------------------------------------------------
# Data structures & constants
# -----------------------------------------------------------------------------
CandidateScore = namedtuple(
    "CandidateScore", "pbs_id ppt_ratio domains canonical score"
)

# BLASTX DB settings
DB_PREFIX    = "./internal_db/ltr_internal_db"
DOMAIN_FASTA = "./data/LTRRT_internal_domains.fasta"
EVALUE_CUT   = 1e-3
MIN_LEN_AA   = 50

# Domain list
DOMAINS = ["GAG", "PR", "RT", "INT", "RH"]

# Scoring weights (total 100)
PBS_W    = 30   # PBS identity ×30
PPT_W    = 10   # PPT purity ×10
DOM_W    = 8    # each domain ×8 (max 5 → 40)
CANON_W  = 20   # bonus if domains in canonical order

# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def hamming(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def load_trna_tails(trna_fasta: str, tail_len: int = 18):
    tails = {}
    for rec in SeqIO.parse(trna_fasta, "fasta"):
        tail = str(rec.seq)[-tail_len:].upper()
        tails[rec.id] = str(BioSeq(tail).reverse_complement())
    return tails

def detect_pbs(window: str, trna_tails, min_len=12, max_len=18, max_mismatches=3):
    best = (False, None, 0.0)
    w = window.upper()
    for trna_id, rc_tail in trna_tails.items():
        for L in range(max_len, min_len - 1, -1):
            pat = rc_tail[-L:]
            for i in range(len(w) - L + 1):
                seg = w[i:i+L]
                mism = hamming(seg, pat)
                if mism <= max_mismatches:
                    ident = (L - mism) / L
                    if ident > best[2]:
                        best = (True, trna_id, ident)
    return best

def is_ppt(seg: str, purine_threshold=0.7):
    s = seg.upper()
    pur = sum(1 for b in s if b in "AG")
    return pur/len(s) >= purine_threshold and "TTT" not in s and "CCC" not in s

def detect_ppt(window: str, min_len=10, max_len=20, purine_threshold=0.7):
    best = 0.0
    w = window.upper()
    for L in range(max_len, min_len - 1, -1):
        for i in range(len(w) - L + 1):
            seg = w[i:i+L]
            if is_ppt(seg, purine_threshold):
                ratio = (seg.count("A") + seg.count("G")) / L
                if ratio > best:
                    best = ratio
    return (best > 0), best

def make_blast_db(fasta: str, db_prefix: str) -> None:
    if not Path(f"{db_prefix}.pin").exists():
        subprocess.check_call([
            "makeblastdb", "-in", fasta,
            "-dbtype", "prot", "-parse_seqids",
            "-out", db_prefix
        ])

def detect_domains_blastx(nuc_seq: str):
    with tempfile.TemporaryDirectory() as tmp:
        qf = Path(tmp) / "qry.fa"
        outf = Path(tmp) / "out.tsv"
        qf.write_text(f">I\n{nuc_seq}\n")
        subprocess.run([
            "blastx", "-query", str(qf), "-db", DB_PREFIX,
            "-evalue", str(EVALUE_CUT),
            "-outfmt", "6 sseqid length pident qstart qend bitscore",
            "-max_target_seqs", "20", "-seg", "yes"
        ], stdout=open(outf, "w"), check=True)

        found = set()
        for line in outf.read_text().splitlines():
            sseqid, length, pident, *_ = line.split()
            if int(length) < MIN_LEN_AA or float(pident) < 25:
                continue
            dom_full = sseqid.split("__", 1)[0]
            dom = dom_full.split("-", 1)[1] if "-" in dom_full else dom_full
            if dom in DOMAINS:
                found.add(dom)
        return [d for d in DOMAINS if d in found]

# -----------------------------------------------------------------------------
# Scoring truncated candidates
# -----------------------------------------------------------------------------
def score_truncated(seq, trna_tails):
    found_pbs, _, pbs_id = detect_pbs(str(seq[:40]), trna_tails)
    pbs_score = (pbs_id if found_pbs else 0.0) * PBS_W

    found_ppt, ppt_ratio = detect_ppt(str(seq[-40:]))
    ppt_score = (ppt_ratio if found_ppt else 0.0) * PPT_W

    doms = detect_domains_blastx(str(seq))
    domain_score = len(doms) * DOM_W

    canon = (doms == [d for d in DOMAINS if d in doms])
    canon_score = CANON_W if canon else 0

    total = pbs_score + ppt_score + domain_score + canon_score
    return CandidateScore(pbs_id, ppt_ratio, doms, canon, total)

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-fasta',  required=True)
    p.add_argument('--candidates',   required=True,
                   help="TSV from LTR pair detector")
    p.add_argument('--trna',         required=True)
    p.add_argument('--output-tsv',   required=True)
    p.add_argument('--output-fasta', required=True)
    p.add_argument('--score-threshold', type=float, default=50.0)
    args = p.parse_args()

    make_blast_db(DOMAIN_FASTA, DB_PREFIX)
    trna_tails = load_trna_tails(args.trna)
    idx = SeqIO.index(args.input_fasta, 'fasta')

    excluded = set()
    with open(args.candidates) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seqid = row['seq_id'].split(":", 1)[0]
            excluded.add(seqid)

    results = []
    count = 0
    for seq_id, rec in idx.items():
        count += 1
        print("processing candidate:" + str(count))
        if seq_id in excluded:
            continue

        cs = score_truncated(rec.seq, trna_tails)
        if cs.score < args.score_threshold:
            continue

        # capture domain coords
        with tempfile.TemporaryDirectory() as tmp:
            qf = Path(tmp) / "qry.fa"
            outf = Path(tmp) / "out_coords.tsv"
            qf.write_text(f">{seq_id}\n{rec.seq}\n")
            subprocess.run([
                "blastx","-query",str(qf),"-db",DB_PREFIX,
                "-evalue",str(EVALUE_CUT),
                "-outfmt","6 qstart qend sseqid",
                "-max_target_seqs","20","-seg","yes"
            ], stdout=open(outf,"w"), check=True)

            starts, ends = [], []
            for line in outf.read_text().splitlines():
                qstart, qend, sseqid = line.split()
                dom_full = sseqid.split("__",1)[0]
                dom = dom_full.split("-",1)[1] if "-" in dom_full else dom_full
                if dom in DOMAINS:
                    starts.append(int(qstart))
                    ends.append(int(qend))

            if not starts:
                continue

            # define range robustly
            positions = starts + ends
            flank = 20
            start = max(0, min(positions) - flank)
            end   = min(len(rec.seq), max(positions) + flank)

        results.append({
            'seq_id':       seq_id,
            'PBS_identity': f"{cs.pbs_id*100:.1f}",
            'PPT_ratio':    f"{cs.ppt_ratio*100:.1f}",
            'domains':      ",".join(cs.domains) or "-",
            'canonical':    int(cs.canonical),
            'trunc_score':  f"{cs.score:.1f}",
            'start':        start,
            'end':          end
        })

    # write TSV
    fields = ['seq_id','PBS_identity','PPT_ratio','domains',
              'canonical','trunc_score','start','end']
    with open(args.output_tsv, 'w', newline='') as out:
        w = csv.DictWriter(out, fields, delimiter='\t')
        w.writeheader()
        for r in results:
            w.writerow(r)

    # write FASTA with headers "<seq_id>":<start>_<end>
    records = []
    for r in results:
        rec = idx[r['seq_id']]
        subseq = rec.seq[r['start']:r['end']]
        header = f"{r['seq_id']}:{r['start']}_{r['end']}"
        records.append(SeqRecord(subseq, id=header, description=""))
    SeqIO.write(records, args.output_fasta, 'fasta')

    print(f"Wrote {len(results)} truncated entries to {args.output_tsv}")
    print(f"Wrote {len(records)} sequences to {args.output_fasta}")

if __name__ == '__main__':
    main()

