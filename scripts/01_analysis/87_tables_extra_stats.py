# 87_tables_extra_stats.py — 给 S2a/S6a/S8 补上"第二套统计"列，并更新 S3c 说明
# 说明：审稿意见指出原表只给单一检验口径（t 检验 / 单侧方向 p / binomial 0.5），
#       需同时给出 donor-level Wilcoxon、双向 BH、以及 Fisher 背景校正，避免"选择性报告"。
import os, sys, io
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from openpyxl import load_workbook

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
TAB = os.path.join(ROOT, "submission", "Tables")
R = lambda *p: os.path.join(ROOT, *p)

NOTE_HDR = "NOTE"


def rewrite_sheet(path, sheet, df, note):
    """用 df 覆盖某个 sheet（保留 NOTE 横幅 + 表头），其余 sheet 不动。"""
    wb = load_workbook(path)
    idx = wb.sheetnames.index(sheet)
    del wb[sheet]
    ws = wb.create_sheet(sheet, idx)
    ws.cell(row=1, column=1, value=NOTE_HDR)
    ws.cell(row=2, column=1, value=note)
    for j, c in enumerate(df.columns, start=1):
        ws.cell(row=4, column=j, value=str(c))
    for i, row in enumerate(df.itertuples(index=False), start=5):
        for j, v in enumerate(row, start=1):
            if isinstance(v, (np.floating, float)):
                if np.isnan(v):
                    continue
                v = float(f"{v:.6g}")
            elif isinstance(v, np.integer):
                v = int(v)
            ws.cell(row=i, column=j, value=v)
    wb.save(path)
    print(f"  [{sheet}] {df.shape[0]} 行 × {df.shape[1]} 列")


# ---------------------------------------------------------------- S2a
p2 = R("celltype_proportion_stats.csv")
d2 = pd.read_csv(p2)
d2 = d2.rename(columns={"padj": "padj_ttest"})
d2 = d2[["celltype", "mean_young", "mean_aged", "fc", "log2FC",
         "p_ttest", "padj_ttest", "p_wilcox", "padj_wilcox"]]
rewrite_sheet(
    os.path.join(TAB, "Table_S2_Cell_composition_statistics.xlsx"),
    "S2a_proportion_tests", d2,
    "Table S2a. Donor-level composition tests (arcsine-sqrt transformed fractions, 5 vs 5). "
    "Both a Welch t-test (p_ttest) and an exact Wilcoxon rank-sum test (p_wilcox) are reported; "
    "padj_ttest and padj_wilcox are Benjamini-Hochberg adjustments across the cell types tested "
    "within each test. The two tests agree in direction and, for SMC/pericyte and EC, in "
    "significance; the FAP increase does not reach significance under either test "
    "(p_ttest = 0.083, p_wilcox = 0.095), which is stated in the text.")

# ---------------------------------------------------------------- S6a
p6 = R("pathway_summary_by_celltype.csv")
d6 = pd.read_csv(p6).rename(columns={"padj2": "padj_dir_two_sided"})
d6 = d6[["celltype", "gene_set", "n_genes_in_set", "n_up", "n_down", "mean_lfc",
         "median_lfc", "frac_up", "p_up", "p_down", "direction", "p_dir",
         "padj_dir", "padj_dir_two_sided"]]
rewrite_sheet(
    os.path.join(TAB, "Table_S6_Pathway_and_clock_genes.xlsx"),
    "S6a_pathway_summary", d6,
    "Table S6a. Directional pathway summary per cell type (15 literature-curated gene sets). "
    "p_dir and padj_dir use one-sided directional tests (BH across 15 sets per cell type), the "
    "manner in which these sets were pre-specified. padj_dir_two_sided additionally charges the "
    "correction for having tested both directions (BH across 30 tests per cell type) and is the "
    "conservative value quoted where both are given; the clock and SASP/ECM results in FAPs "
    "survive this stricter correction.")

# ---------------------------------------------------------------- S8
p8 = R("regulatory_chain_FAP.csv")
d8 = pd.read_csv(p8)
N_BG, P0 = 16096, 0.493           # 细胞类型背景：16,096 个可测基因中 49.3% 为正 log2FC
bg_up, bg_dn = int(round(N_BG * P0)), N_BG - int(round(N_BG * P0))


def fisher_row(r):
    a = int(round(r.frac_rna_up * r.n_genes_measured))
    b = int(r.n_genes_measured) - a
    orr, p = stats.fisher_exact([[a, b], [bg_up, bg_dn]], alternative="two-sided")
    return pd.Series({"n_up_measured": a, "n_down_measured": b,
                      "background_frac_up": P0, "fisher_OR": orr, "p_fisher": p})


d8 = pd.concat([d8, d8.apply(fisher_row, axis=1)], axis=1)
d8["padj_fisher"] = multipletests(d8["p_fisher"], method="fdr_bh")[1]
d8 = d8[["TF", "DA_direction", "n_DA_peaks_with_motif", "n_linked_genes",
         "n_genes_measured", "n_up_measured", "n_down_measured", "frac_rna_up",
         "background_frac_up", "median_rna_lfc", "p_binom", "fisher_OR", "p_fisher",
         "padj_fisher", "top_genes"]]
d8.to_csv(p8, index=False)
rewrite_sheet(
    os.path.join(TAB, "Table_S8_Regulatory_chain_FAP.xlsx"),
    "S8_regulatory_chain", d8,
    "Table S8. TF -> enhancer -> target-gene regulatory chains in aged FAPs. Per TF and direction: "
    "number of motif-bearing DA peaks, number of linked genes measured, fraction up-regulated, the "
    "exact binomial p against 0.5 (p_binom), and a Fisher exact test of the 2x2 table "
    "(motif-linked targets vs all other expressed genes) x (up vs down) using the FAP cell-type "
    "background of 49.3% positive log2FC among 16,096 expressed genes (fisher_OR, p_fisher, "
    "padj_fisher; BH across all 48 TF-direction rows). Because the background is not 50%, the "
    "Fisher test is the appropriate primary test; p_binom is retained for transparency.")

# ---------------------------------------------------------------- S3c 说明
p3 = os.path.join(TAB, "Table_S3c_All_DA_DESeq2_results.xlsx")
wb = load_workbook(p3, read_only=True)
ws = wb[wb.sheetnames[0]]
old_note = ws.cell(row=2, column=1).value
wb.close()
new_note = (
    "Table S3c. All DESeq2 differential-accessibility results (5 vs 5, pseudo-bulk per donor; "
    "45,713 support-10 consensus regions x 6 cell types, 228,987 rows). Two effect-size thresholds "
    "are used in the paper: (i) padj < 0.05 alone, giving 1,588 DA peaks in total (FAP 1,122 "
    "[821 open, 73%]; myonuclei 331; MuSC 82; EC 49; SMC/pericyte 3; macrophage 1); and (ii) "
    "padj < 0.05 with |log2FC| > 1, the threshold used for all downstream motif, coupling and "
    "regulatory-chain analyses and plotted in Fig. 2A, giving 836 DA peaks (FAP 408 [310 open, "
    "76%]; myonuclei 318; MuSC 70; EC 37; SMC/pericyte 2; macrophage 1). Fig. 2A bars show the "
    "second set; peak numbers quoted for FAPs in the Abstract and Results refer to the first set "
    "where explicitly stated.")
wb = load_workbook(p3)
wb[wb.sheetnames[0]].cell(row=2, column=1, value=new_note)
wb.save(p3)
print("  [S3c] 说明已更新（两套阈值计数写入表注）")
print("\nS8 前 3 行预览：")
print(d8.head(3)[["TF", "DA_direction", "n_genes_measured", "frac_rna_up",
                  "p_binom", "p_fisher", "padj_fisher"]].to_string(index=False))
