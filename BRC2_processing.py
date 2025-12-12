from Bio.Seq import Seq
from collections import defaultdict
import ahocorasick

# Reference sequence
REF = ("ATGTTCTTGCAGAGGAGAACAAAAATAGTGTAAAGCAGCATATAAAAATGACTCTAGGTCAAGATTTAAAATCGGACATCTCCTTGAATATAGATAAAATACCAGAAAAAAATAATGATTACATGAACAAATGGGCAGGACTCTTAGGTCCAATTTCAAATGAAGTGGGATTTAGGGGCTTTTATTCTGCTCATGGCACAAAACTGAATGTTTCTACTGAAGCTCTGCAAAAAGCTGTGAAACTGTTTAGTGATATTGAGAATATCAGAACTGAGCATAGTCTTCACTATTCACCTACGTCTAGACAAAATGTATCAAAAATACTTCCTCGTGTTGATAAGAGAAACCCAGAGCACTGTGTAAACT")

regions = [179, 182, 185, 188]  # 0-based
FASTQ = "XBRC2Sat-DMSORP2.merged.fastq"
OUTFILE = "variant_counts_BRSat_DMSO_RP2.tsv"

# Create all possible 3-mer mutations for given regions
mutants = {}
for pos in regions:
    wt = REF[pos:pos+3]
    for b1 in "ACGT":
        for b2 in "ACGT":
            for b3 in "ACGT":
                alt = b1 + b2 + b3
                if alt == wt:
                    continue
                mut_seq = REF[:pos] + alt + REF[pos+3:]
                rc_seq = str(Seq(mut_seq).reverse_complement())
                label = f"{pos+1}-X-{alt}"  # 1-based
                mutants[label] = (mut_seq, rc_seq)

# Build Aho-Corasick automaton for fast multiple pattern matching
automaton = ahocorasick.Automaton()
for label, (mut_seq, rc_seq) in mutants.items():
    automaton.add_word(mut_seq, (label, 'F'))
    automaton.add_word(rc_seq, (label, 'R'))
automaton.make_automaton()

# Initialize counts
counts = defaultdict(int)

# Process FASTQ file
with open(FASTQ) as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            seq = line.strip()
            for end_idx, (label, _) in automaton.iter(seq):
                counts[label] += 1

# Write output
with open(OUTFILE, "w") as out:
    out.write("Mutation\tCount\n")
    for label in sorted(mutants):
        out.write(f"{label}\t{counts[label]}\n")

print(f"Done. Output written to {OUTFILE}")
