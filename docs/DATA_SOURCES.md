# Data sources

All data analysed here are publicly available; no new human material was collected.
Download each resource, then run the scripts in the order given in `ANALYSIS_ORDER.md`.

| Resource | What it is | Accession / location | Used by |
|---|---|---|---|
| GSE268953 | Single-nucleus multiome (RNA + ATAC) of human hamstring muscle, 5 young and 5 aged male donors | NCBI GEO, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268953 | 01–18 |
| GTEx v10 | Bulk RNA-seq of skeletal muscle, 818 donors | https://gtexportal.org | 15, 18 |
| GSE132040 | Tabula Muris Senis bulk RNA-seq, mouse limb muscle across age | GEO | 20, 21 |
| PXD043715 | Mouse gastrocnemius proteome, ageing time course | PRIDE | 22 |
| PXD015616 | Human muscle proteome (Ubaida-Mohien et al., eLife 2019), 58 donors | PRIDE; per-channel tables in the paper's supplementary files | 24, 26 |
| OMIX004308 / OMIX004305 | Lai et al., Nature 2024 limb-muscle atlas: processed object and per-cell-type snATAC peak matrices | National Genomics Data Center (NGDC) OMIX | 67, 68, 73, 74 |
| — | Kedlian et al. intercostal-muscle ageing atlas (processed objects) | https://www.muscleageingcellatlas.org; raw data ArrayExpress E-MTAB-13874 | 66 |
| — | Lehallier et al., Nat Med 2019 plasma proteome ageing (SomaScan, 4,263 people) | supplementary tables of the paper | 59–61 |
| — | Human Protein Atlas plasma concentration | https://www.proteinatlas.org | 60 |
| hg38 refGene | Transcript models used for peak annotation and locus plots | UCSC `refGene` table via https://api.genome.ucsc.edu | 11, 16, `fetch_gene_coords.py` |

## Not redistributed here

Raw sequencing data, intermediate Seurat/ArchR objects and the processed result
tables are **not** included in this repository (size). Processed tables are provided
as Supplementary Tables S1–S14 of the manuscript.
