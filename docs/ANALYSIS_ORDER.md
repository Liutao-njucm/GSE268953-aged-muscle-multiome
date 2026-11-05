# Analysis order

Stages are sequential; numbers within a stage are the script prefixes.
`R` scripts are run with `Rscript`, Python scripts with `python`.

| Stage | Scripts | What happens | Main outputs |
|---|---|---|---|
| 1. Nucleus-level object | `01`–`04b` | Load 10 samples, split GEX/ATAC, scDblFinder doublets, Harmony integration, clustering, marker-based annotation | annotated multiome object; cell-type proportions |
| 2. Chromatin | `05*`–`09b` | Per-sample chromVAR, consensus peak set, DA peaks at donor level (support = 10), motif enrichment | DA peak tables, motif enrichment (Table S3, S4) |
| 3. Transcription | `10`–`13b` | Donor-level pseudobulk DESeq2 per cell type, pathway summary, clock-gene panel | DEG tables (Table S5), pathway table (Table S6) |
| 4. Coupling and regulation | `11`–`14` | Peak-to-gene assignment, chromatin–transcriptome coupling classes, TF → enhancer → target chains | coupling tables (Table S7), regulatory chain (Table S8) |
| 5. External validation | `15`, `17`–`18c`, `20`–`24a` | GTEx linear models, module-wise concordance, mouse bulk and proteome comparisons | Tables S9–S11 |
| 6. Regulon activity | `29a`–`30b` | CollecTRI + decoupler ULM, threshold-free target-set shift | Table S12 |
| 7. Translatability and replication | `54`, `59`–`61`, `65`–`68`, `73`–`74` | Plasma proteome projection, Lai 2024 snRNA/snATAC replication, detection-rate gate | Tables S13, S14 |
| 8. Figures | `16`, `19`, `21`, `22`, `26`, `31`, `45`, `47`–`49`, `55`, `57`–`58`, `61`, `68b`, `80`, `104`, `112` | Panel rendering at 300 and 600 dpi, final composition | Fig. 1–5, Fig. S1–S5 |

Scripts that assemble the manuscript text, render the docx and run textual
consistency checks are **not** included: they operate on local intermediate files
and do not contribute to any result or figure. They are available from the
corresponding author on request.
