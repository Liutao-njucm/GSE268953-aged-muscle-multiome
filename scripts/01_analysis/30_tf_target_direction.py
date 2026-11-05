"""
30_tf_target_direction.py

Directional test of TF target sets against cell-type DEG tables.
Complements the (underpowered) donor-level ULM regulon scores in 29_*: here each
TF's CollecTRI target set is tested against the DEG direction calls, which uses
per-gene DESeq2 evidence and is therefore far better powered at n = 5 vs 5.

For every TF with >= 10 measured targets:
  * restrict to targets that are significant DEGs (padj < 0.05, |LFC| > 1)
  * binomial test of the up-regulated fraction against the cell-type background
    (fraction of all significant DEGs that are up)
  * Fisher exact test on the 2x2 table (target vs non-target) x (up vs down)

Output: results/GSE268953/regulon/tf_target_direction.csv
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

CELLTYPES = ["FAP", "Myonuclei", "MuSC", "EC", "T_NK_cells", "Macrophage", "SMC_Pericyte"]

net = pd.read_csv(NET, sep="\t")
net = net[["source_genesymbol", "target_genesymbol"]].dropna().drop_duplicates()
net.columns = ["TF", "target"]
net["TF"] = net["TF"].str.upper()
net["target"] = net["target"].str.upper()
net = net[net["TF"] != net["target"]]

rows = []
for ct in CELLTYPES:
    f = os.path.join(DEG, f"DEG_{ct}.csv")
    if not os.path.exists(f):
        print("missing", f); continue
    d = pd.read_csv(f)
    d["gene"] = d["gene"].astype(str).str.upper()
    d = d.drop_duplicates("gene")
    sig = d[(d["padj"] < 0.05) & (d["log2FoldChange"].abs() > 1)].copy()
    sig["dir"] = np.where(sig["log2FoldChange"] > 0, "up", "down")
    p0 = (sig["dir"] == "up").mean()
    up_genes = set(sig.loc[sig["dir"] == "up", "gene"])
    dn_genes = set(sig.loc[sig["dir"] == "down", "gene"])
    n_up_tot, n_dn_tot = len(up_genes), len(dn_genes)

    for tf, sub in net.groupby("TF"):
        tg = set(sub["target"])
        if len(tg & set(d["gene"])) < 10:
            continue
        a = len(tg & up_genes)
        b = len(tg & dn_genes)
        if a + b < 5:
            continue
        # background: significant DEGs not in this target set
        rest_up = n_up_tot - a
        rest_dn = n_dn_tot - b
        if rest_up + rest_dn < 10:
            continue
        orr, pf = stats.fisher_exact([[a, b], [rest_up, rest_dn]], alternative="two-sided")
        pb = stats.binomtest(a, a + b, p0, alternative="two-sided").pvalue
        rows.append({"celltype": ct, "TF": tf, "n_targets_measured": int(len(tg & set(d["gene"]))),
                     "n_sig_targets": a + b, "n_up": a, "n_down": b,
                     "frac_up": a / (a + b), "background_frac_up": p0,
                     "binom_p": pb, "fisher_OR": orr, "fisher_p": pf})

res = pd.DataFrame(rows)
res["padj"] = np.nan
for ct, idx in res.groupby("celltype").groups.items():
    res.loc[idx, "padj"] = multipletests(res.loc[idx, "fisher_p"], method="fdr_bh")[1]
res = res.sort_values(["celltype", "fisher_p"])
res.to_csv(os.path.join(REG, "tf_target_direction.csv"), index=False)

CLOCK = ["NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "HLF", "TEF", "ARNTL", "ARNTL2",
         "CLOCK", "NPAS2", "BHLHE40", "BHLHE41", "CIART", "NFIL3"]
AP1 = ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "BATF", "BATF3", "MAFK"]
panel = set(CLOCK + AP1 + ["CEBPA", "CEBPB", "CEBPD", "NFIA", "TEAD1", "MEF2C", "PPARG", "RUNX2"])

print(f"cell types: {res['celltype'].nunique()}   TF-celltype tests: {len(res)}")
print("\n=== FAP: key TFs ===")
k = res[(res["celltype"] == "FAP") & (res["TF"].isin(panel))].sort_values("fisher_p")
print(k.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
print("\n=== FAP: top 15 any TF ===")
print(res[res["celltype"] == "FAP"].head(15).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
print("\n=== clock TFs, all cell types, padj<0.05 ===")
c = res[(res["TF"].isin(CLOCK)) & (res["padj"] < 0.05)]
print(c.sort_values(["celltype", "padj"]).to_string(index=False, float_format=lambda v: f"{v:.3g}") if len(c) else "none")
print("\nDONE")
