import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert LTR Detector output to LTRharvest .scn format"
    )
    p.add_argument("-i", "--input",  required=True,
                   help="LTR Detector .bed‑style file")
    p.add_argument("-o", "--output", required=True,
                   help="Output .scn file")
    return p.parse_args()

def convert_from_tail(seqid, tail, seq_nr):
    """
    tail: list of the last 17 whitespace‐separated fields:
      [0] s_ret, [1] e_ret,
      [2] s_l,   [3] e_l,
      [4] s_r,   [5] e_r,
      [6] sim%,  [7] TG_s, [8] TG_e,
      [9] CA_s,  [10] CA_e,
      [11] PPT_s,[12] PPT_e,
      [13] strand,[14] purine%,
      [15] TSD_s,[16] TSD_e
    """
    # full element
    s_ret, e_ret = int(tail[0]), int(tail[1])
    l_ret = e_ret - s_ret + 1

    # left LTR
    s_l, e_l = int(tail[2]), int(tail[3])
    l_l = e_l - s_l + 1

    # right LTR
    s_r, e_r = int(tail[4]), int(tail[5])
    l_r = e_r - s_r + 1

    # always parse sim as float
    sim = float(tail[6].rstrip('%'))

    return [
        s_ret, e_ret, l_ret,
        s_l,   e_l,   l_l,
        s_r,   e_r,   l_r,
        sim,       
        seq_nr,
        seqid
    ]

def main():
    args = parse_args()
    seq_map = {}
    next_seq_nr = 0

    with open(args.input) as fin, open(args.output, "w") as fout:
        fout.write("# Converted from LtrDetector\n")
        fout.write("# s(ret) e(ret) l(ret) "
                   "s(lLTR) e(lLTR) l(lLTR) "
                   "s(rLTR) e(rLTR) l(rLTR) "
                   "sim seq-nr chr\n")

        for line in fin:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 18:
                continue

            tail = parts[-17:]
            if not tail[0].isdigit():
                continue

            seqid = parts[0]
            if seqid not in seq_map:
                seq_map[seqid] = next_seq_nr
                next_seq_nr += 1

            rec = convert_from_tail(seqid, tail, seq_map[seqid])

            rec[9] = f"{rec[9]:.1f}"
            fout.write(" ".join(str(x) for x in rec) + "\n")

if __name__ == "__main__":
    main()
