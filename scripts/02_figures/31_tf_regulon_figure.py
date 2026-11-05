"""
31_tf_regulon_figure.py

Figure S8: transcriptional read-out of TF activity (4 panels).
A  heatmap of target-set log2FC shift (median target - median rest) across cell types
B  FAP volcano of target-set shift vs -log10 p
C  FAP target-gene log2FC distributions for representative TFs vs background
D  donor-level ULM regulon scores (honest null result at n = 5 vs 5)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = "D:/workbuddy/outputs/2026-10-25-20-23-31"
REG = os.path.join(ROOT, "results/GSE268953/regulon")
DEG = os.path.join(ROOT, "results/GSE268953/deg")
FIG = os.path.join(ROOT, "results/GSE268953/figures")
os.makedirs(FIG, exist_ok=True)

sh = pd.read_csv(os.path.join(REG, "tf_target_lfc_shift.csv"))
ulm = pd.read_csv(os.path.join(REG, "regulon_activity_by_donor.csv"))
dif = pd.read_csv(os.path.join(REG, "regulon_activity_differential.csv"))

CLOCK = ["NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "HLF", "TEF", "ARNTL", "NPAS2",
         "CLOCK", "BHLHE40"]
AP1 = ["JUN", "JUNB", "JUND", "FOS", "FOSL1", "FOSL2", "BATF", "MAFK"]
OTH = ["CEBPA", "CEBPB", "CEBPD", "NFIA", "TEAD1", "MEF2C", "PPARG", "RUNX2"]
PANEL = CLOCK + AP1 + OTH
CTS = ["FAP", "Myonuclei", "MuSC", "EC"]

fig = plt.figure(figsize=(14.5, 10.5))
gs = fig.add_gridspec(2, 2, wspace=0.42, hspace=0.42,
                      left=0.13, right=0.97, top=0.93, bottom=0.08)

# ---------------------------------------------------------------- A heatmap
axA = fig.add_subplot(gs[0, 0])
mat = np.full((len(PANEL), len(CTS)), np.nan)
pl = np.ones((len(PANEL), len(CTS)))
for i, tf in enumerate(PANEL):
    for j, ct in enumerate(CTS):
        r = sh[(sh.TF == tf) & (sh.celltype == ct)]
        if len(r):
            mat[i, j] = r.delta_median_lfc.values[0]
            pl[i, j] = r.p.values[0]
vmax = 0.30
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im = axA.imshow(np.clip(mat, -vmax, vmax), cmap="RdBu_r", norm=norm, aspect="auto")
axA.set_xticks(range(len(CTS)), [c.replace("_", "/") for c in CTS], fontsize=10)
axA.set_yticks(range(len(PANEL)), PANEL, fontsize=8.5)
for i in range(len(PANEL)):
    for j in range(len(CTS)):
        if np.isfinite(mat[i, j]) and pl[i, j] < 0.05:
            axA.text(j, i, "*", ha="center", va="center", fontsize=13, color="k")
axA.axhline(len(CLOCK) - 0.5, color="k", lw=0.8, ls="--")
axA.axhline(len(CLOCK) + len(AP1) - 0.5, color="k", lw=0.8, ls="--")
axA.set_title("a  TF target-set shift\n(median log2FC of targets - rest)", fontsize=11, loc="left")
cb = fig.colorbar(im, ax=axA, fraction=0.03, pad=0.02)
cb.set_label("Δ median log2FC", fontsize=9); cb.ax.tick_params(labelsize=8)

# ---------------------------------------------------------------- B volcano (FAP)
axB = fig.add_subplot(gs[0, 1])
f = sh[sh.celltype == "FAP"].copy()
f["nlp"] = -np.log10(f.p)
sig = f.p < 0.05
axB.scatter(f.loc[~sig, "delta_median_lfc"], f.loc[~sig, "nlp"], s=12, c="#B4B2A9", alpha=0.6, lw=0)
axB.scatter(f.loc[sig, "delta_median_lfc"], f.loc[sig, "nlp"], s=26, c="#D85A30", alpha=0.9, lw=0)
lab = f[(f.p < 0.01) | (f.TF.isin(["NR1D1", "CLOCK", "DBP"]))].copy()
lab = lab.sort_values(["nlp", "delta_median_lfc"])
xlim = axB.get_xlim()
placed = []
for _, r in lab.iterrows():
    x0, y0 = r.delta_median_lfc, r.nlp
    dy = 0.0
    for px, py in placed:
        if abs(px - x0) < 0.12 and abs(py - (y0 + dy)) < 0.16:
            dy += 0.17
    placed.append((x0, y0 + dy))
    ## 靠近右边缘的标签放点的左侧，避免超出画布被截断
    if x0 > 0.55 * xlim[1]:
        axB.annotate(r.TF, (x0, y0), fontsize=7.5, ha="right",
                     xytext=(x0 - 0.02, y0 + dy), textcoords="data")
    else:
        axB.annotate(r.TF, (x0, y0), fontsize=7.5, xytext=(x0 + 0.02, y0 + dy),
                     textcoords="data")
axB.axhline(-np.log10(0.05), color="k", lw=0.7, ls=":")
axB.axvline(0, color="k", lw=0.6)
axB.set_xlabel("Δ median log2FC (targets - rest)", fontsize=10)
axB.set_ylabel("-log10 P (Mann-Whitney)", fontsize=10)
axB.set_title("b  FAP target-set directionality\n(481 TFs tested)", fontsize=11, loc="left")
## FDR 注释放左上空白区，避开右上高 P 散点标签
axB.text(0.02, 0.98, "transcriptome-wide FDR: none < 0.05\nPPARG padj$_{panel}$ = 0.026 (26-TF panel)\nnominal p shown",
         transform=axB.transAxes, ha="left", va="top", fontsize=8.5, color="#5F5E5A")

# ---------------------------------------------------------------- C distributions
axC = fig.add_subplot(gs[1, 0])
net = pd.read_csv(os.path.join(ROOT, "data/regulon/collectri_omnipath_2026.tsv"), sep="\t")
net = net[["source_genesymbol", "target_genesymbol"]].dropna().drop_duplicates()
net.columns = ["TF", "target"]
net["TF"] = net.TF.str.upper(); net["target"] = net.target.str.upper()
d = pd.read_csv(os.path.join(DEG, "DEG_FAP.csv"))
d["gene"] = d.gene.astype(str).str.upper()
d = d.dropna(subset=["log2FoldChange"]).drop_duplicates("gene")
lfc = d.set_index("gene")["log2FoldChange"]

show_tf = ["NR1D1", "CLOCK", "FOSL2", "CEBPB", "PPARG"]
data, labels, pv = [lfc.values], ["rest"], [np.nan]
for tf in show_tf:
    tg = [g for g in net.loc[net.TF == tf, "target"] if g in lfc.index]
    data.append(lfc[lfc.index.isin(tg)].values)
    p = sh[(sh.TF == tf) & (sh.celltype == "FAP")].p
    labels.append(f"{tf}\n(n={len(tg)})" + (f"\np={p.values[0]:.3f}" if len(p) else ""))
    pv.append(p.values[0] if len(p) else np.nan)
bp = axC.boxplot(data, vert=False, widths=0.6, showfliers=False, patch_artist=True,
                 tick_labels=labels, medianprops=dict(color="k", lw=1.2))
## 面板色统一：所有 TF 箱线同色（rest 灰）——原紫/橙区分无生物学含义且图注未解释
for patch, c in zip(bp["boxes"], ["#D3D1C7"] + ["#F0997B"] * len(show_tf)):
    patch.set_facecolor(c); patch.set_alpha(0.8)
axC.axvline(0, color="k", lw=0.7, ls=":")
axC.set_xlabel("log2FC (aged / young), FAP", fontsize=10)
axC.tick_params(axis="y", labelsize=8)
axC.set_title("c  FAP target-gene log2FC distributions", fontsize=11, loc="left")

# ---------------------------------------------------------------- D ULM per donor
axD = fig.add_subplot(gs[1, 1])
show = ["NR1D1", "CLOCK", "JUNB", "FOSL2"]
pos = np.arange(len(show))
for k, tf in enumerate(show):
    sub = ulm[(ulm.TF == tf) & (ulm.celltype == "FAP")]
    for gi, (grp, col, dx) in enumerate([("Young", "#378ADD", -0.16), ("Aged", "#D85A30", 0.16)]):
        v = sub[sub.group == grp].score.values
        axD.scatter(np.full(len(v), k + dx) + np.random.default_rng(1).normal(0, 0.02, len(v)),
                    v, s=30, c=col, alpha=0.9, lw=0, zorder=3)
        axD.hlines(np.median(v), k + dx - 0.07, k + dx + 0.07, color="k", lw=1.4, zorder=4)
    pj = dif[(dif.TF == tf) & (dif.celltype == "FAP")].padj
    if len(pj):
        ## padj 标签放在面板内顶部（数据坐标），与散点最高点(3.2)分离，也不挤标题
        axD.annotate(f"padj={pj.values[0]:.2f}", (k, 3.5),
                     ha="center", fontsize=8, color="#5F5E5A", zorder=5)
axD.set_xticks(pos, show, fontsize=10)
axD.set_ylim(-5.2, 3.9)
axD.set_ylabel("ULM regulon score", fontsize=10)
axD.set_title("d  Donor-level regulon activity (FAP)\nno TF survives FDR at n = 5 vs 5", fontsize=11, loc="left")
axD.legend(handles=[Line2D([], [], marker="o", ls="", color="#378ADD", label="Young"),
                    Line2D([], [], marker="o", ls="", color="#D85A30", label="Aged")],
           fontsize=8.5, frameon=False, loc="lower right")
axD.axhline(0, color="k", lw=0.6, ls=":")

fig.suptitle("Transcriptional read-out of TF activity: target-set directionality is detectable, "
             "donor-level regulon scores are not", fontsize=12, x=0.02, ha="left", y=0.985)
out = os.path.join(FIG, "FigS8_TF_regulon_activity.png")
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print("saved", out)
