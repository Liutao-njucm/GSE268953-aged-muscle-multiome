# Script inventory

Auto-generated from the header comment of each script.

### 01_analysis — data processing and statistics

Run in the order given in `ANALYSIS_ORDER.md`.

| Script | Purpose |
|---|---|
| `01_load_qc.R` | GSE268953 snMultiome: load 10 samples, split GEX/ATAC, per-sample QC |
| `02_doublet.R` | per-sample doublet removal via scDblFinder (Seurat5-compatible) |
| `03_integrate.R` | merge, Harmony integration, clustering, marker genes |
| `04_annotate.R` | rule-based scoring annotation + group proportions + plots |
| `04b_fix_labels.R` | 04b_fix_labels.R |
| `05_chromvar.R` | Phase 2A - Step 1: chromVAR TF motif activity on ATAC |
| `05a_rebuild_atac.R` | Phase 2A - Step 0 (v2): rebuild multiome object, fast path |
| `05b_chromvar_per_sample.R` | Phase 2A - Step 1 (v3): per-sample chromVAR, peak-capped, incremental |
| `05c_consensus_da.R` | Phase 2A - Step 2: consensus-peak donor-level DA (Aged vs Young) |
| `05d_merge_chromvar.R` | Phase 2A - Step 1d: merge per-sample chromVAR results |
| `05f_atac_pipeline.R` | 05f: ATAC per-sample processing (chromVAR + consensus pseudobulk) |
| `06_aging_clock.R` | Phase 2A - Step 3: chromatin aging clock |
| `06_merge_chromvar.R` | 06: merge per-sample chromVAR deviations into one matrix |
| `07_chromatin_aging_clock.R` | 07: chromatin aging clock analysis |
| `08_consensus_DA.R` | 08: differential chromatin accessibility (consensus peaks) |
| `08a_consensus_support.R` | 08a: consensus peak sample-support |
| `09_motif_enrichment.R` | 09: motif enrichment on differential accessibility peaks |
| `09b_map_motifs.R` | 09b: map JASPAR matrix IDs -> TF names, rebuild aging signature |
| `10_rna_pseudobulk_de.R` | 10: RNA donor-level pseudobulk differential expression (5 vs 5) |
| `11_peak_annotation.R` | 11: annotate consensus peaks to genes (hg38, UCSC refGene) |
| `12_coupling_analysis.R` | 12: chromatin-transcriptome coupling (per cell type) |
| `13_pathway_summary.R` | 13: pathway-level summary of the ageing transcriptome |
| `13b_clock_genes.R` |  |
| `14_regulatory_chain.R` | 14: TF -> DA peak -> target gene regulatory chain (FAP focus) |
| `15_gtex_validation.py` | OLS with intercept + one binary covariate; returns (b_group, p_group) via t-test. |
| `17_proportion_stats.R` | 17: cell-type proportion test (5v5 donor-level, arcsine-sqrt) |
| `18_concordance.R` | 18: cross-dataset concordance — snRNA FAP log2FC vs GTEx bulk |
| `18b_fap_clock_donorplot.R` | 18b: per-donor FAP expression of clock genes (robustness plot) |
| `18c_concordance_modules.R` | 18c: module-wise GTEx concordance (honest split: SASP/ECM vs clock vs other) |
| `20_tms_clock_check.py` | Tabula Muris Senis bulk RNA-seq, limb muscle across age |
| `24_human_muscle_proteome.py` | 24: donor-level analysis of the Ubaida-Mohien et al. eLife 2019 human skeletal muscle TMT proteome (58 donors, 20-34 -> 80+ yr), and concordance with  |
| `24a_map_uniprot.py` | 24a: map UniProt entry names (e.g. 1433B_HUMAN) -> gene symbols. Uses UniProt REST batch queries (id: query) because the stream endpoint times out her |
| `29_tf_regulon_activity.py` | 29_tf_regulon_activity.py TF regulon activity from donor-level pseudobulk RNA, using CollecTRI (OmniPath) and the Univariate Linear Model (ULM; decoup |
| `29a_prep_regulon_input.R` | 29a_prep_regulon_input.R |
| `30_tf_target_direction.py` | 30_tf_target_direction.py Directional test of TF target sets against cell-type DEG tables. Complements the (underpowered) donor-level ULM regulon scor |
| `30b_tf_target_lfc_shift.py` | 30b_tf_target_lfc_shift.py Threshold-free complement to 30_*: for each TF, compare the log2FC distribution of its CollecTRI target genes against all n |
| `38_build_supp_tables.py` | Build submission-ready supplementary tables Table S1-S12 as multi-sheet xlsx. Each workbook is named by its in-text identifier so figure/table callout |
| `46_split_table_s3.py` | 46_split_table_s3.py Split the 44 MB Table_S3_DA_peaks_and_annotation.xlsx into three standalone workbooks (S3a / S3b / S3c) because BMC caps individu |
| `54_add_prg4.R` | compute per-donor FAP CPM10k for PRG4 (+ validation genes) |
| `59_plasma_biomarker.py` | map the FAP aging signature onto a population plasma proteome. Plasma source : Lehallier et al., Nat Med 2019 (SomaScan, 4,263 donors, 2,925 assays; l |
| `59a_explore_plasma.py` | explore intersection of the FAP SASP/ECM signature with the Lehallier 2019 plasma proteome (SomaScan, 4,263 donors, 2,925 proteins). |
| `59b_plasma_gsea.py` | classic weighted Kolmogorov-Smirnov running sum (unweighted hits). |
| `60_plasma_translatability.py` | build the plasma translatability table for the FAP SASP/ECM program. For every gene in our SASP / ECM / clock / muscle modules: * FAP snRNA log2FC / p |
| `65_pilot_gse167186.R` | GSE167186 pilot: detection-rate gate for clock repressor arm in FAP |
| `66_pilot_kedlian.py` | Kedlian Nat Aging 2024 stromal/endothelial compartments |
| `67_pilot_lai.py` | Lai et al. Nature 2024 (OMIX004308-02 stromal h5ad, limb muscle 15-99y) |
| `73_inspect_lai_atac.R` | 解包 OMIX004305-05 并检查 per-celltype peak matrix RDS 结构 |
| `74_lai_atac_clock.R` | Lai OMIX004305 snATAC 时钟位点可及性复制分析 |
| `87_tables_extra_stats.py` | 用 df 覆盖某个 sheet（保留 NOTE 横幅 + 表头），其余 sheet 不动。 |
| `fetch_gene_coords.py` |  |
| `install_annotation.R` |  |
| `install_deseq2.R` |  |

### 02_figures — figure rendering

Each script writes one panel or one composed figure.

| Script | Purpose |
|---|---|
| `104_regchain_only.R` | 只重绘 Fig 4D（regulatory chain）面板 |
| `112_final_polish.R` | 全图终审后的三个面板重绘（2026-11-05） |
| `16_locus_plots.R` | 16: locus accessibility tracks for clock genes |
| `19_main_figures.R` | Phase 3: auxiliary panels + assemble 5 main figures |
| `21_tms_clock_figure.py` | TMS mouse limb-muscle clock figure |
| `22_pxd043715_concordance.py` | 碰撞感知标签：依次尝试候选偏移，选第一个与已放置标签无重叠的位置； 靠近右缘的点优先把文字放到点左侧，避免跨入右子图。 |
| `26_fig5_proteome_panel.py` | Fig5 panel B: human muscle proteome — ECM/SASP protein module (Ubaida-Mohien eLife 2019). Left: ECM protein log2FC (old vs young) with FDR stars; Righ |
| `31_tf_regulon_figure.py` | 31_tf_regulon_figure.py Figure S8: transcriptional read-out of TF activity (4 panels). A heatmap of target-set log2FC shift (median target - median re |
| `45_export_600dpi.R` | re-render all Fig.1-5 panels at 600 dpi |
| `47_export_600dpi_compose.py` | 47_export_600dpi_compose.py 600-dpi re-export of the five main figures, exactly mirroring the layout of scripts/27_assemble_main_figures.py (auto-trim |
| `47a_fig5_proteome_panel_600dpi.py` | 47a_fig5_proteome_panel_600dpi.py Re-render the Fig5 proteome panel (scripts/26) at 600 dpi into panels_600dpi/. Does NOT overwrite figures/fig5_prote |
| `48_fix_panels.R` | redraw the five problem panels |
| `49_compose_final.py` | assemble the five main figures from panels_600dpi/ at two scales: scale 1 (3000 px @ 300 dpi) -> submission/Figures/Figure1..5.png (submission) + figu |
| `55_redraw_5c.R` | 55 (v2): redraw Fig 5C clock-locus tracks, genome-browser style |
| `57_teaching_nr1d1.R` | 57 (v2): teaching figure for Fig 5C — NR1D1 locus, numbered callouts |
| `58_redraw_5a.R` | redraw only the Fig5A GTEx concordance panel |
| `61_plot_plasma_fig.py` | Supplementary Figure S9: plasma translatability of the FAP SASP/ECM program. (A) evidence dot-matrix: FAP log2FC / GTEx log2FC / plasma age coefficien |
| `68_figS5_lai_replication.py` | independent replication of clock repressor arm in Lai Nat 2024 |
| `68b_figS5_600dpi.py` | 68_figS5_lai_replication.py --- independent replication of clock repressor arm in Lai Nat 2024 |
| `80_fix_fig1c.R` | 重绘 Fig. 1C 面板：修正副标题中错误的 p 值 |
