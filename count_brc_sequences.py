#!/usr/bin/env python3
"""Count reads containing each full-length BRC sequence (forward and reverse complement)."""

import argparse
import gzip

from Bio import SeqIO
from Bio.Seq import Seq


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fastq", help="merged FASTQ (.fastq or .fastq.gz)")
    p.add_argument("-s", "--seqs", default="refs/BRC_all.fa", help="FASTA of target sequences")
    p.add_argument("-o", "--out", help="output TSV (default: stdout)")
    a = p.parse_args()

    targets = {r.id: str(r.seq).upper() for r in SeqIO.parse(a.seqs, "fasta")}
    rc = {k: str(Seq(v).reverse_complement()) for k, v in targets.items()}
    fwd = dict.fromkeys(targets, 0)
    rev = dict.fromkeys(targets, 0)
    total = 0

    opener = gzip.open if a.fastq.endswith(".gz") else open
    with opener(a.fastq, "rt") as f:
        for rec in SeqIO.parse(f, "fastq"):
            total += 1
            read = str(rec.seq)
            for k in targets:
                fwd[k] += targets[k] in read
                rev[k] += rc[k] in read

    lines = ["ID\tForward\tReverse\tTotal"]
    lines += [f"{k}\t{fwd[k]}\t{rev[k]}\t{fwd[k] + rev[k]}" for k in targets]
    lines.append(f"# total_reads\t{total}")
    text = "\n".join(lines) + "\n"

    if a.out:
        open(a.out, "w").write(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
