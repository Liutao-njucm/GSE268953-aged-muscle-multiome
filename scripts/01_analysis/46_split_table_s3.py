"""
46_split_table_s3.py
Split the 44 MB Table_S3 workbook into three parts (S3a/S3b/S3c) for smaller uploads.
"""
# Rebuilds each part from the source CSVs (identical content to script 38),
# leaving the original combined workbook untouched. The two largest parts are
# additionally provided as compressed CSV.
"""
import os
import pandas as pd

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
OUT = os.path.join(ROOT, "submission", "Tables")
os.makedirs(OUT, exist_ok=True)
R = lambda *p: os.path.join(ROOT, *p)

print("reading source CSVs ...")
da = pd.read_csv(R("DA_peaks_Aged_vs_Young_by_celltype.csv"))
ann = pd.read_csv(R("peak_annotation.csv"))
sig = da[(da["padj"] < 0.05) & (da["log2FoldChange"].abs() > 1)].copy()
ann_slim = ann[["peak_id", "chr", "start", "end", "width", "support", "symbol",
                "dist_to_tss", "in_gene_body", "region_class", "linked_gene"]].copy()
sig = sig.merge(ann_slim[["peak_id", "symbol", "dist_to_tss", "region_class", "linked_gene"]],
                on="peak_id", how="left")
sig = sig.sort_values(["celltype", "padj"])
print(f"  full DA rows = {len(da):,}; significant = {len(sig):,}; peaks annotated = {len(ann_slim):,}")


def write_book(fname, sheets):
    path = os.path.join(OUT, fname)
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for title, obj, note in sheets:
            if note:
                pd.DataFrame({"NOTE": [note]}).to_excel(
                    xw, sheet_name=title[:31], index=False, startrow=0)
                obj.to_excel(xw, sheet_name=title[:31], index=False, startrow=3)
            else:
                obj.to_excel(xw, sheet_name=title[:31], index=False)
    mb = os.path.getsize(path) / 1e6
    print(f"  {fname}  {mb:6.1f} MB")
    return mb


# ---- S3a: significant DA peaks -------------------------------------------
write_book("Table_S3a_Significant_DA_peaks.xlsx", [
    ("S3a_significant_DA_peaks", sig,
     "Table S3a. Significant differentially accessible peaks (padj<0.05, |log2FC|>1; "
     "donor pseudo-bulk DESeq2, 5 vs 5) with gene annotation; supports Fig. 2A. "
     "Full unfiltered DESeq2 output: Table S3c."),
])

# ---- S3b: annotation of all consensus peaks -------------------------------
write_book("Table_S3b_Consensus_peak_annotation.xlsx", [
    ("S3b_peak_annotation", ann_slim,
     "Table S3b. Annotation of all 283,971 consensus peaks "
     "(UCSC hg38 refGene; TSS within 100 kb or gene body)."),
])

# ---- S3c: complete DESeq2 output ------------------------------------------
mb_c = write_book("Table_S3c_All_DA_DESeq2_results.xlsx", [
    ("S3c_all_DA_results", da,
     "Table S3c. Complete DESeq2 output for all tested peaks, all cell types."),
])

# ---- compressed CSV fallback for S3c (in case xlsx exceeds upload cap) ----
csv_gz = os.path.join(OUT, "Table_S3c_All_DA_DESeq2_results.csv.gz")
da.to_csv(csv_gz, index=False, compression="gzip")
print(f"  Table_S3c_All_DA_DESeq2_results.csv.gz  {os.path.getsize(csv_gz)/1e6:6.1f} MB")

print("\nDone. Original combined workbook kept unchanged:")
print("  ", os.path.join(OUT, "Table_S3_DA_peaks_and_annotation.xlsx"),
      f"{os.path.getsize(os.path.join(OUT, 'Table_S3_DA_peaks_and_annotation.xlsx'))/1e6:.1f} MB")
