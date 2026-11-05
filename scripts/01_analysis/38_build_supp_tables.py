"""
Build submission-ready supplementary tables Table S1-S12 as multi-sheet xlsx.
Each workbook is named by its in-text identifier so figure/table callouts match files.
"""
import os
import pandas as pd

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
OUT = os.path.join(ROOT, "submission", "Tables")
os.makedirs(OUT, exist_ok=True)

R = lambda *p: os.path.join(ROOT, *p)


def sheet_name(s):
    s = s.replace("/", "_").replace("\\", "_")[:31]
    return s


def write_book(fname, sheets):
    """sheets: list of (sheet_title, DataFrame or path-to-csv, note)"""
    path = os.path.join(OUT, fname)
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        first = True
        for title, obj, note in sheets:
            df = pd.read_csv(obj) if isinstance(obj, str) else obj
            if note:
                # a one-row provenance banner above the header block
                banner = pd.DataFrame({"NOTE": [note]})
                banner.to_excel(xw, sheet_name=sheet_name(title), index=False, startrow=0)
                df.to_excel(xw, sheet_name=sheet_name(title), index=False, startrow=3)
            else:
                df.to_excel(xw, sheet_name=sheet_name(title), index=False)
            first = False
    mb = os.path.getsize(path) / 1e6
    print(f"  {fname}  {mb:6.1f} MB")
    return path


# ---------------------------------------------------------------- S1
write_book("Table_S1_Sample_QC_and_cell_composition.xlsx", [
    ("S1a_sample_QC", R("per_sample_qc.csv"),
     "Table S1a. Per-sample QC after doublet removal (snRNA + snATAC)."),
    ("S1b_nuclei_per_sample", R("celltype_counts.csv"),
     "Table S1b. Nuclei retained per donor and cell type."),
    ("S1c_fraction_per_sample", R("celltype_proportions.csv"),
     "Table S1c. Cell-type fraction per donor."),
])

# ---------------------------------------------------------------- S2
write_book("Table_S2_Cell_composition_statistics.xlsx", [
    ("S2a_proportion_tests", R("celltype_proportion_stats.csv"),
     "Table S2a. Donor-level composition tests (arcsine-sqrt transformed fractions, 5 vs 5)."),
    ("S2b_chromatin_aging_score", R("chromatin_aging_score_by_celltype.csv"),
     "Table S2b. ATAC-based chromatin ageing score per cell type (see Fig. 2C)."),
])

# ---------------------------------------------------------------- S3
print("  [S3] reading large peak tables ...")
da = pd.read_csv(R("DA_peaks_Aged_vs_Young_by_celltype.csv"))
ann = pd.read_csv(R("peak_annotation.csv"))
sig = da[(da["padj"] < 0.05) & (da["log2FoldChange"].abs() > 1)].copy()
ann_slim = ann[["peak_id", "chr", "start", "end", "width", "support", "symbol",
                "dist_to_tss", "in_gene_body", "region_class", "linked_gene"]].copy()
sig = sig.merge(ann_slim[["peak_id", "symbol", "dist_to_tss", "region_class", "linked_gene"]],
                on="peak_id", how="left")
sig = sig.sort_values(["celltype", "padj"])
print(f"      full DA rows = {len(da):,}; significant = {len(sig):,}; peaks annotated = {len(ann_slim):,}")
write_book("Table_S3_DA_peaks_and_annotation.xlsx", [
    ("S3a_significant_DA_peaks", sig,
     "Table S3a. Significant differentially accessible peaks (padj<0.05, |log2FC|>1; donor pseudo-bulk DESeq2, 5 vs 5) with gene annotation. Full unfiltered DESeq2 output: Table S3c."),
    ("S3b_peak_annotation", ann_slim,
     "Table S3b. Annotation of all 283,971 consensus peaks (UCSC hg38 refGene; TSS within 100 kb or gene body)."),
    ("S3c_all_DA_results", da,
     "Table S3c. Complete DESeq2 output for all tested peaks, all cell types."),
])

# ---------------------------------------------------------------- S4
write_book("Table_S4_Motif_enrichment_and_TF_activity.xlsx", [
    ("S4a_motif_enrichment", R("motif_enrichment_DA_peaks.csv"),
     "Table S4a. JASPAR2024 motif enrichment in DA peaks by cell type and direction (Fisher, BH)."),
    ("S4b_chromVAR_by_celltype", R("motif_aging_by_celltype_5v5.csv"),
     "Table S4b. chromVAR motif activity, aged vs young, per cell type (5 vs 5)."),
    ("S4c_chromVAR_global", R("motif_aging_global_5v5.csv"),
     "Table S4c. chromVAR motif activity, global (all nuclei)."),
    ("S4d_aging_TF_signature", R("aging_TF_signature.csv"),
     "Table S4d. The 80-TF ageing signature used in Fig. 2B."),
])

# ---------------------------------------------------------------- S5
write_book("Table_S5_DEG_per_celltype.xlsx", [
    ("S5a_DEG_summary", R("deg", "DEG_summary_by_celltype.csv"),
     "Table S5a. DEG counts per cell type (padj<0.05, |log2FC|>1)."),
    ("S5b_all_DEGs", R("deg", "DEG_all_celltypes.csv"),
     "Table S5b. Complete DESeq2 output for all genes and cell types."),
])

# ---------------------------------------------------------------- S6
write_book("Table_S6_Pathway_and_clock_genes.xlsx", [
    ("S6a_pathway_summary", R("pathway_summary_by_celltype.csv"),
     "Table S6a. Directional hypergeometric pathway summary (15 curated gene sets)."),
    ("S6b_pathway_gene_lfc", R("pathway_gene_lfc_long.csv"),
     "Table S6b. Per-gene log2FC within each pathway and cell type."),
    ("S6c_clock_genes", R("clock_genes_by_celltype.csv"),
     "Table S6c. Circadian genes per cell type (drives Fig. 3C)."),
    ("S6d_clock_per_donor", R("FAP_clock_per_donor.csv"),
     "Table S6d. Per-donor FAP expression of clock and SASP genes (5 vs 5 donor separation)."),
])

# ---------------------------------------------------------------- S7
write_book("Table_S7_Chromatin_RNA_coupling.xlsx", [
    ("S7a_coupling_gene_level", R("coupling_gene_level_by_celltype.csv"),
     "Table S7a. Gene-level chromatin-RNA coupling (linked peak score vs RNA log2FC)."),
    ("S7b_coupling_summary", R("coupling_summary_by_celltype.csv"),
     "Table S7b. Coupling class fractions and Fisher tests per cell type (Fig. 4B)."),
    ("S7c_clock_locus_ATAC", R("clock_locus_accessibility.csv"),
     "Table S7c. Pseudo-bulk accessibility at clock gene loci (Fig. 5C)."),
])

# ---------------------------------------------------------------- S8
write_book("Table_S8_Regulatory_chain_FAP.xlsx", [
    ("S8_regulatory_chain", R("regulatory_chain_FAP.csv"),
     "Table S8. TF -> open enhancer -> target-gene chains in aged FAPs (exact binomial test; Fig. 4D)."),
])

# ---------------------------------------------------------------- S9
write_book("Table_S9_GTEx_validation.xlsx", [
    ("S9a_GTEx_candidate_genes", R("gtex", "gtex_muscle_validation.csv"),
     "Table S9a. GTEx v10 skeletal muscle (818 donors) aged vs young log2TPM for candidate genes (Fig. 5A)."),
    ("S9b_GTEx_samples", R("gtex", "gtex_muscle_samples.csv"),
     "Table S9b. The 818 GTEx skeletal-muscle samples analysed (age, sex, batch)."),
    ("S9c_module_concordance", R("gtex_concordance_by_module.csv"),
     "Table S9c. Directional concordance between GTEx bulk and FAP snRNA by module."),
    ("S9d_snRNA_concordance", R("gtex_snRNA_concordance.csv"),
     "Table S9d. Gene-level GTEx versus cell-type-resolved snRNA concordance."),
])

# ---------------------------------------------------------------- S10
write_book("Table_S10_Cross_species_TMS_PXD043715.xlsx", [
    ("S10a_TMS_clock", R("tms_mouse_muscle_clock.csv"),
     "Table S10a. Tabula Muris Senis bulk limb muscle (n=54, 1-27 months), age effect per circadian gene (Fig. S1A)."),
    ("S10b_PXD043715", R("pxd043715_protein_concordance.csv"),
     "Table S10b. PXD043715 mouse gastrocnemius proteome versus human snRNA log2FC (Fig. S1B)."),
])

# ---------------------------------------------------------------- S11
write_book("Table_S11_Human_muscle_proteome.xlsx", [
    ("S11a_proteome_aging", R("proteome", "human_muscle_proteome_aging.csv"),
     "Table S11a. Human vastus lateralis TMT proteome, old (65+) vs young (20-34); re-analysis of Ubaida-Mohien et al. eLife 2019, 58 donors, 5,891 proteins (Fig. 5B, S2)."),
    ("S11b_panel_genes", R("proteome", "proteome_panel_genes.csv"),
     "Table S11b. Curated ECM / clock / myofibre protein panels."),
    ("S11c_concordance_by_celltype", R("proteome", "proteome_concordance_by_celltype.csv"),
     "Table S11c. Bulk-proteome versus cell-type-resolved snRNA concordance (Fig. S2F)."),
    ("S11d_module_scores", R("proteome", "proteome_module_scores_by_donor.csv"),
     "Table S11d. Per-donor module scores (ECM, clock, myofibre)."),
])

# ---------------------------------------------------------------- S12
write_book("Table_S12_TF_regulon_activity.xlsx", [
    ("S12a_target_lfc_shift", R("regulon", "tf_target_lfc_shift.csv"),
     "Table S12a. CollecTRI target-set log2FC shift per TF and cell type (Mann-Whitney); padj = transcriptome-wide BH, padj_panel = pre-specified 26-TF panel (Fig. S3B-C)."),
    ("S12b_target_direction", R("regulon", "tf_target_direction.csv"),
     "Table S12b. Directional Fisher test of target sets among significant DEGs."),
    ("S12c_ULM_differential", R("regulon", "regulon_activity_differential.csv"),
     "Table S12c. Donor-level ULM regulon activity, aged vs young (Fig. S3D)."),
    ("S12d_panel_key_TFs", R("regulon", "regulon_panel_key_tfs.csv"),
     "Table S12d. The pre-specified clock / AP-1 / CEBP / other hypothesis panel."),
])

print("\nDone. Tables written to", OUT)
