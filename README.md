# BRC repeats saturation editing

Read counting and heatmap for BRC repeats saturation variants.

## Scripts

| Script | What it does |
|---|---|
| `scripts/count_saturation_variants.py` | Counts reads for every codon substitution at the target positions for BRC2 |
| `scripts/count_brc_sequences.py` | Counts reads matching each full BRC sequence |
| `scripts/BRC2_heatmap.R` | Plots log2(ratio) heatmap per amino acid and treatment |

Reference sequences are in `refs/`.

## Usage

```bash
pip install biopython pyahocorasick

python scripts/count_saturation_variants.py sample.merged.fastq -o sample_counts.tsv
python scripts/count_brc_sequences.py sample.merged.fastq -o sample_brc.tsv
Rscript scripts/BRC2_heatmap.R counts_BRC2.xlsx BRC2_heatmap.pdf
```

R packages: `openxlsx`, `dplyr`, `tidyr`, `ComplexHeatmap`, `circlize`.
