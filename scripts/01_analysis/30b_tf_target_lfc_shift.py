"""
30b_tf_target_lfc_shift.py

Threshold-free complement to 30_*: for each TF, compare the log2FC distribution of
its CollecTRI target genes against all non-target genes (Mann-Whitney U). This keeps
small regulons usable (e.g. NR1D1 has ~12 targets) and does not depend on an
arbitrary DEG cutoff, so it is the best-powered transcriptional read-out available
from 5 vs 5 donors.

Output: results/GSE268953/regulon/tf_target_lfc_shift.csv
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = "D:/workbuddy/outputs/2026-10-25-20-23-31"
DEG = os.path.join(ROOT, "results/GSE268953/deg")
REG = os.path.join(ROOT, "results/GSE268953/regulon")
NET = os.path.join(ROOT, "data/regulon/collectri_omnipath_2026.tsv")

CELLTYPES = ["FAP", "Myonuclei", "MuSC", "EC"]

net = pd.read_csv(NET, sep="\t")[["source_genesymbol", "target_genesymbol"]].dropna().drop_duplicates()
net.columns = ["TF", "target"]
net["TF"] = net["TF"].str.upper()
net["target"] = net["target"].str.upper()
net = net[net["TF"] != net["target"]]

rows = []
for ct in CELLTYPES:
    d = pd.read_csv(os.path.join(DEG, f"DEG_{ct}.csv"))
    d["gene"] = d["gene"].astype(str).str.upper()
    d = d.dropna(subset=["log2FoldChange"]).drop_duplicates("gene")
    lfc = d.set_index("gene")["log2FoldChange"]
    for tf, sub in net.groupby("TF"):
        tg = [g for g in sub["target"] if g in lfc.index]
        if len(tg) < 10:
            continue
        in_set = lfc.index.isin(tg)
        a = lfc[in_set].values
        b = lfc[~in_set].values
        if len(b) < 100:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"celltype": ct, "TF": tf, "n_targets": len(tg),
                     "median_lfc_target": float(np.median(a)),
                     "median_lfc_rest": float(np.median(b)),
                     "delta_median_lfc": float(np.median(a) - np.median(b)),
                     "mw_u": u, "p": p})

res = pd.DataFrame(rows)
res["padj"] = np.nan
for ct, idx in res.groupby("celltype").groups.items():
    res.loc[idx, "padj"] = multipletests(res.loc[idx, "p"], method="fdr_bh")[1]
res = res.sort_values(["celltype", "p"])

# ---- pre-specified hypothesis panel (clock / AP-1 / CEBP / adipogenic / myogenic TFs).
# Correction within this panel is reported alongside the transcriptome-wide correction,
# because these TFs were specified a priori from the biological model, not selected
# post hoc from the top of the list.
CLOCK = ["NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "HLF", "TEF", "ARNTL", "ARNTL2",
         "CLOCK", "NPAS2", "BHLHE40", "BHLHE41", "CIART", "NFIL3"]
AP1 = ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "BATF", "BATF3", "MAFK"]
OTHER = ["CEBPA", "CEBPB", "CEBPD", "NFIA", "TEAD1", "MEF2C", "PPARG", "RUNX2", "SMAD3", "STAT3"]
PANEL = set(CLOCK + AP1 + OTHER)

res["padj_panel"] = np.nan
for ct, idx in res[res.TF.isin(PANEL)].groupby("celltype").groups.items():
    res.loc[idx, "padj_panel"] = multipletests(res.loc[idx, "p"], method="fdr_bh")[1]
res["n_panel_tested"] = res.groupby("celltype")["padj_panel"].transform(lambda s: s.notna().sum())
res.to_csv(os.path.join(REG, "tf_target_lfc_shift.csv"), index=False)

CLOCK = ["NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "HLF", "TEF", "ARNTL", "ARNTL2",
         "CLOCK", "NPAS2", "BHLHE40", "BHLHE41", "CIART", "NFIL3"]
AP1 = ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "BATF", "BATF3", "MAFK"]
OTHER = ["CEBPA", "CEBPB", "CEBPD", "NFIA", "TEAD1", "MEF2C", "PPARG", "RUNX2", "SMAD3", "STAT3"]

for ct in CELLTYPES:
    s = res[res["celltype"] == ct]
    k = s[s["TF"].isin(CLOCK + AP1 + OTHER)].sort_values("p")
    print(f"\n=== {ct}: key TFs, {len(s)} tested ===")
    print(k.head(12).to_string(index=False, float_format=lambda v: f"{v:.3g}"))

print("\n=== clock TFs with padj<0.05 (all cell types) ===")
c = res[(res["TF"].isin(CLOCK)) & (res["padj"] < 0.05)]
print(c.sort_values(["celltype", "padj"]).to_string(index=False, float_format=lambda v: f"{v:.3g}") if len(c) else "none")
print("\nDONE")
