# Analysis code — Circadian repressor collapse primes the aged muscle niche

Analysis code for the manuscript of the same title; citation details will be added once the
article is published.

The study integrates single-nucleus ATAC and RNA sequencing of human hamstring muscle
from five young and five aged male donors (GSE268953) to ask which regulatory elements
change with age inside each cell type, which transcription factors rewrite them, and how
much of that chromatin remodelling reaches the transcriptome.

**This repository contains code only.** Raw data are public (see `docs/DATA_SOURCES.md`);
processed result tables are provided as Supplementary Tables S1–S14 of the manuscript.

## Layout

```
scripts/01_analysis/     data processing, statistics and result tables
scripts/02_figures/      figure rendering (300 / 600 dpi)
docs/DATA_SOURCES.md     every external dataset, with accession
docs/ANALYSIS_ORDER.md   run order, stage by stage
docs/SCRIPT_INVENTORY.md one-line description of every script
env/                     R sessionInfo and Python package list
```

## Getting started

1. Install R 4.5.2 and Python 3.13; the exact package versions used are in `env/`.
2. Download the resources listed in `docs/DATA_SOURCES.md` and place them under `data/`.
3. Run every script **from the repository root**: paths inside the scripts are relative
   (`data/...` for input, `results/GSE268953/...` for output).
4. Run the stages in the order given in `docs/ANALYSIS_ORDER.md`.

Scripts are numbered by stage. Numbers are preserved from the working project, so they are
not contiguous; `docs/ANALYSIS_ORDER.md` is the authoritative order.

## Statistical conventions

* Peaks and genes are called at **donor level**, never at nucleus level; the donor is the
  unit of replication throughout (n = 5 per group).
* Differential accessibility uses only the 45,713 consensus regions supported by all ten
  donors (`support = 10`), which removes the missing-peak bias of per-sample peak sets.
* Two-sided tests with BH correction; where a pre-specified panel was used, both the
  panel-wide and the transcriptome-wide correction are reported.
* Effect sizes and confidence intervals accompany every test statistic quoted in the paper.

## Licence

Code: MIT (`LICENSE`). Processed tables in the manuscript: CC BY 4.0.

## Contact

Tao Liu, School of Chinese Medicine, Nanjing University of Chinese Medicine —
njzyywhx@126.com
