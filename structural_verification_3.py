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
Candidate = namedtuple(
    "Candidate",
    "seq_id ltr1_start ltr1_end ltr2_start ltr2_end identity score gap support"
)


# BLASTX-based domain detection config
DB_PREFIX    = "./ltr_internal_db/ltr_internal_db"
DOMAIN_FASTA = "./data/LTRRT_internal_domains.fasta"
EVALUE_CUT   = 1e-3
MIN_LEN_AA   = 50  # minimum aligned AA length for domain

# Domain scoring
DOMAINS = ["GAG","PR","RT","INT","RH"]
DOMAIN_POINTS = 5         # per domain
CANON_BONUS   = 5         # if domains appear in canonical order

# -----------------------------------------------------------------------------
# PBS detection
# -----------------------------------------------------------------------------
def load_trna_tails(trna_fasta: str, tail_len: int = 18):
    """Load tRNA 3' tails and store their reverse complements."""
    tails = {}
    for rec in SeqIO.parse(trna_fasta, "fasta"):
        tail = str(rec.seq)[-tail_len:].upper()
        tails[rec.id] = str(BioSeq(tail).reverse_complement())
    return tails

def hamming(a: str, b: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(a, b))

def detect_pbs(window: str, trna_tails, min_len=12, max_len=18, max_mismatches=3):
    """
    Scan `window` for tRNA-tail matches.
    Returns (found, tRNA_id, identity_frac [0–1]).
    """
    best = (False, None, 0.0)
    w = window.upper()
    for trna_id, rc_tail in trna_tails.items():
        for L in range(max_len, min_len -1, -1):
            pat = rc_tail[-L:]
            for i in range(len(w) - L + 1):
                seg = w[i:i+L]
                mism = hamming(seg, pat)
                if mism <= max_mismatches:
                    ident = (L - mism) / L
                    if ident > best[2]:
                        best = (True, trna_id, ident)
    return best

# -----------------------------------------------------------------------------
# PPT detection
# -----------------------------------------------------------------------------
def is_ppt(seg: str, purine_threshold=0.7):
    s = seg.upper()
    pur = sum(1 for b in s if b in "AG")
    return pur/len(s) >= purine_threshold and "TTT" not in s and "CCC" not in s

def detect_ppt(window: str, min_len=10, max_len=20, purine_threshold=0.7):
    """
    Scan `window` for the highest-purity purine-rich tract.
    Returns (found, purity_ratio [0–1]).
    """
    best = 0.0
    w = window.upper()
    for L in range(max_len, min_len -1, -1):
        for i in range(len(w) - L + 1):
            seg = w[i:i+L]
            if is_ppt(seg, purine_threshold):
                ratio = (seg.count("A") + seg.count("G")) / L
                if ratio > best:
                    best = ratio
    return (best > 0), best

# -----------------------------------------------------------------------------
# BLASTX-based internal-domain detection
# -----------------------------------------------------------------------------
def make_blast_db(fasta: str, db_prefix: str) -> None:
    """Build (or skip) a protein BLAST DB from `fasta`."""
    if not Path(f"{db_prefix}.pin").exists():
        subprocess.check_call([
            "makeblastdb", "-in", fasta,
            "-dbtype", "prot", "-parse_seqids",
            "-out", db_prefix
        ])

def detect_domains_blastx(nuc_seq: str) -> list:
    """
    BLASTX the `nuc_seq` against DB_PREFIX and return
    the detected domains (GAG, PR, RT, INT, RH) in canonical order.
    """
    with tempfile.TemporaryDirectory() as tmp:
        qf  = Path(tmp) / "qry.fa"
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
            dom_full = sseqid.split("__",1)[0]    # e.g. "Ty3-RT"
            dom = dom_full.split("-",1)[1] if "-" in dom_full else dom_full
            if dom in DOMAINS:
                found.add(dom)

        return [d for d in DOMAINS if d in found]

def compute_features(seq, five_start, five_end, three_start, identity, trna_tails,
                     support, max_support, support_points):
    
    # 1) LTR similarity → 0–30
    sim_score = max(0.0, min(identity,1.0)) * 30.0

    # 2) PBS → 0–20
    w1 = str(seq[five_end:five_end+40])
    pbs_found, pbs_trna, pbs_id = detect_pbs(w1, trna_tails)
    pbs_score = (pbs_id if pbs_found else 0.0) * 20.0

    # 3) PPT → 0–10
    w2 = str(seq[max(0, three_start-40):three_start])
    ppt_found, ppt_ratio = detect_ppt(w2)
    ppt_score = ppt_ratio * 10.0

    # 4) Domains → 5 points each
    internal = str(seq[five_end:three_start])
    doms = []
    if len(internal) >= 150:
        doms = detect_domains_blastx(internal)
    domain_score = sum(DOMAIN_POINTS for d in doms)

    # 5) Canonical order → 5
    canon_bonus = CANON_BONUS if doms == [d for d in DOMAINS if d in doms] else 0

    # 6) Support → up to support_points
    #    Scale support linearly to max_support
    sup_frac = min(support, max_support) / max_support
    support_score = sup_frac * support_points

    total = sim_score + pbs_score + ppt_score + domain_score + canon_bonus + support_score

    return {
        'sim_score':     sim_score,
        'pbs_score':     pbs_score,
        'PBS_hit':       pbs_trna or '-',
        'PBS_identity':  pbs_id if pbs_found else 0.0,
        'ppt_score':     ppt_score,
        'PPT_hit':       str(ppt_found),
        'PPT_ratio':     ppt_ratio,
        'domains':       doms,
        'domain_score':  domain_score,
        'canon_bonus':   canon_bonus,
        'support_score': support_score,
        'total_score':   total
    }

def analyse_candidate(rec, cand, trna_tails, max_support, support_points):
    # evaluate both orientations
    featA = compute_features(rec.seq,
                             cand.ltr1_start, cand.ltr1_end,
                             cand.ltr2_start, cand.identity,
                             trna_tails, cand.support,
                             max_support, support_points)
    featB = compute_features(rec.seq,
                             cand.ltr2_start, cand.ltr2_end,
                             cand.ltr1_start, cand.identity,
                             trna_tails, cand.support,
                             max_support, support_points)
    chosen = featB if featB['total_score'] > featA['total_score'] else featA

    return {
        'seq_id':        cand.seq_id,
        'ltr1_start':    cand.ltr1_start,
        'ltr1_end':      cand.ltr1_end,
        'ltr2_start':    cand.ltr2_start,
        'ltr2_end':      cand.ltr2_end,
        'ltr_similarity':f"{cand.identity*100:.2f}",
        'PBS_hit':       chosen['PBS_hit'],
        'PBS_identity':  f"{chosen['PBS_identity']*100:.2f}",
        'PPT_hit':       chosen['PPT_hit'],
        'PPT_ratio':     f"{chosen['PPT_ratio']*100:.2f}",
        'domains':       ','.join(chosen['domains']) if chosen['domains'] else '-',
        'support':       cand.support,
        'support_score': f"{chosen['support_score']:.2f}",
        'total_score':   f"{chosen['total_score']:.2f}"
    }

# -----------------------------------------------------------------------------
# I/O and main
# -----------------------------------------------------------------------------
def load_candidates(path):
    c = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for r in reader:
            c.append(Candidate(
                r['seq_id'],
                int(r['ltr1_start']), int(r['ltr1_end']),
                int(r['ltr2_start']), int(r['ltr2_end']),
                float(r['identity']),
                float(r['score']),
                int(r['gap']),
                int(r.get('support', 0))   # new support column
            ))
    return c

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-fasta',  required=True)
    p.add_argument('--candidates',   required=True)
    p.add_argument('--trna',         required=True)
    p.add_argument('--output-tsv',   required=True)
    p.add_argument('--output-fasta', required=True)
    p.add_argument('--score-threshold', type=float, default=70.0)
    # new support parameters:
    p.add_argument('--max-support',   type=int,   default=100,
                   help="Support value that maps to full support-points")
    p.add_argument('--support-points', type=float, default=10.0,
                   help="Maximum points to award for support")
    args = p.parse_args()

    make_blast_db(DOMAIN_FASTA, DB_PREFIX)
    trna_tails = load_trna_tails(args.trna)
    idx = SeqIO.index(args.input_fasta, 'fasta')
    cands = load_candidates(args.candidates)

    results = []
    count = 0
    for cand in cands:
        count += 1
        print("processing candidate:" + str(count))
        if cand.seq_id not in idx:
            continue
        rec = idx[cand.seq_id]
        results.append(analyse_candidate(
            rec, cand, trna_tails,
            args.max_support, args.support_points
        ))

    # sort by total_score desc
    results.sort(key=lambda x: float(x['total_score']), reverse=True)

    # write TSV with new support columns
    fields = [
        'seq_id','ltr1_start','ltr1_end','ltr2_start','ltr2_end',
        'ltr_similarity','PBS_hit','PBS_identity',
        'PPT_hit','PPT_ratio','domains',
        'support','support_score','total_score'
    ]
    with open(args.output_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fields, delimiter='\t')
        writer.writeheader()
        for r in results:
            writer.writerow({k: r[k] for k in fields})

    # write FASTA above threshold
    records = []
    thr = args.score_threshold
    for r in results:
        if float(r['total_score']) < thr:
            break
        rec = idx[r['seq_id']]
        s1,e1 = int(r['ltr1_start']), int(r['ltr1_end'])
        s2,e2 = int(r['ltr2_start']), int(r['ltr2_end'])
        start, end = min(s1,s2), max(e1,e2)
        subseq = rec.seq[start:end]
        header = f"{r['seq_id']}:{start}-{end}"
        records.append(SeqRecord(subseq, id=header, description=""))
    SeqIO.write(records, args.output_fasta, 'fasta')

    print(f"Wrote {len(results)} entries to {args.output_tsv}")
    print(f"Wrote {len(records)} sequences to {args.output_fasta} (score ≥ {thr})")

if __name__ == '__main__':
    main()
