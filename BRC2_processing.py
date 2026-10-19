#!/usr/bin/env python3
"""Count reads carrying each codon substitution at selected positions of the BRC2 reference."""

import argparse
import gzip
from collections import defaultdict
from itertools import product

import ahocorasick
from Bio import SeqIO
from Bio.Seq import Seq


def read_ref(path):
    return str(next(SeqIO.parse(path, "fasta")).seq).upper()


def build_mutants(ref, positions):
    """All 63 non-WT codons at each 0-based position -> {label: full mutant sequence}."""
    mutants = {}
    for pos in positions:
        wt = ref[pos:pos + 3]
        for alt in map("".join, product("ACGT", repeat=3)):
            if alt != wt:
                mutants[f"{pos + 1}-X-{alt}"] = ref[:pos] + alt + ref[pos + 3:]
    return mutants


def build_automaton(mutants):
    auto = ahocorasick.Automaton()
    for label, seq in mutants.items():
        auto.add_word(seq, label)
        auto.add_word(str(Seq(seq).reverse_complement()), label)
    auto.make_automaton()
    return auto


def count(fastq, auto):
    counts = defaultdict(int)
    opener = gzip.open if fastq.endswith(".gz") else open
    with opener(fastq, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                for _, label in auto.iter(line.strip()):
                    counts[label] += 1
    return counts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fastq", help="merged FASTQ (.fastq or .fastq.gz)")
    p.add_argument("-o", "--out", required=True, help="output TSV")
    p.add_argument("-r", "--ref", default="refs/BRC2_ref.fa")
    p.add_argument("-p", "--positions", type=int, nargs="+", default=[179, 182, 185, 188],
                   help="0-based codon start positions")
    a = p.parse_args()

    mutants = build_mutants(read_ref(a.ref), a.positions)
    counts = count(a.fastq, build_automaton(mutants))

    with open(a.out, "w") as out:
        out.write("Mutation\tCount\n")
        for label in sorted(mutants):
            out.write(f"{label}\t{counts[label]}\n")
    print(f"Written: {a.out}")


if __name__ == "__main__":
    main()
