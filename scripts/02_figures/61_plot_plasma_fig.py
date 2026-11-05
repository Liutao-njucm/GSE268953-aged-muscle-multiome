# -*- coding: utf-8 -*-
"""61_plot_plasma_fig.py - Supplementary Figure S9: plasma translatability of the
FAP SASP/ECM program.

(A) evidence dot-matrix: FAP log2FC | GTEx log2FC | plasma age coefficient | plasma q
(B) gene-set enrichment of FAP programs along the plasma age axis (NES)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle, Patch

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
OUT = os.path.join(RES, "plasma")
FIG = os.path.join(RES, "figures")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                     "axes.linewidth": 0.7, "pdf.fonttype": 42})

d = pd.read_csv(os.path.join(OUT, "plasma_translatability.csv"))
g = pd.read_csv(os.path.join(OUT, "plasma_gsea_by_geneset.csv"))
gt = pd.read_csv(os.path.join(RES, "gtex", "gtex_muscle_validation.csv")).rename(
    columns={"p_linear_model": "gtex_p"})
d = d.merge(gt[["gene", "gtex_p"]], on="gene", how="left")

COL = {"SASP_inflammatory": "#E67E22", "ECM_fibrosis": "#27AE60",
       "Clock_repressive_arm": "#8E44AD", "Muscle_fiber": "#2E86C1"}

# ---------------------------------------------------------------- panel A
sub = d[d.module.isin(["SASP_inflammatory", "ECM_fibrosis"])].copy()
# order: module, then FAP log2FC descending
sub["mod_rank"] = sub.module.map({"SASP_inflammatory": 0, "ECM_fibrosis": 1})
sub = sub.sort_values(["mod_rank", "fap_lfc"], ascending=[True, False]).reset_index(drop=True)

cols = [("fap_lfc", "FAP\nlog2FC", 1.6, "fap_padj", 0.05),
        ("gtex_lfc", "GTEx\nlog2FC", 1.2, "gtex_p", 0.05),
        ("plasma_coef_per_decade", "plasma age\ncoef (dec)", 1.0, None, None),
        ("plasma_q", "plasma\nq", 0.0, None, None)]

fig = plt.figure(figsize=(7.4, 11.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1], wspace=0.42,
                      left=0.20, right=0.965, top=0.935, bottom=0.065)
axA = fig.add_subplot(gs[0, 0])
fig.text(0.005, 0.985, "(a) Translatability of the FAP SASP/ECM program to blood",
         fontsize=9.6, fontweight="bold", va="top")

nrow = len(sub)
xpos = [0, 1, 2, 3]
for ci, (key, lab, vmax, padj_col, padj_thr) in enumerate(cols):
    vals = sub[key].values.astype(float)
    if key == "plasma_q":
        # colour by -log10 q, capped
        vv = -np.log10(np.clip(vals, 1e-45, 1))
        norm = TwoSlopeNorm(vmin=0, vcenter=1.3, vmax=np.nanmax(vv))
        cmap = plt.get_cmap("YlOrRd")
        for ri, v in enumerate(vv):
            if np.isnan(v):
                axA.add_patch(Rectangle((ci - .42, ri - .42), .84, .84,
                                        facecolor="#EFEFEF", edgecolor="white", lw=.8))
            else:
                axA.add_patch(Rectangle((ci - .42, ri - .42), .84, .84,
                                        facecolor=cmap(norm(v)), edgecolor="white", lw=.8))
    else:
        lim = np.nanmax(np.abs(vals)) if np.nanmax(np.abs(vals)) > 0 else 1
        cmap = plt.get_cmap("RdBu_r")
        norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
        for ri, v in enumerate(vals):
            if np.isnan(v):
                axA.add_patch(Rectangle((ci - .42, ri - .42), .84, .84,
                                        facecolor="#EFEFEF", edgecolor="white", lw=.8))
            else:
                axA.add_patch(Rectangle((ci - .42, ri - .42), .84, .84,
                                        facecolor=cmap(norm(v)), edgecolor="white", lw=.8))
                if padj_col is not None and pd.notna(sub[padj_col].iloc[ri]):
                    if sub[padj_col].iloc[ri] < padj_thr:
                        axA.add_patch(Rectangle((ci - .44, ri - .44), .88, .88,
                                                facecolor="none", edgecolor="black", lw=1.0))
    axA.text(ci, -1.4, lab, ha="center", va="bottom", fontsize=7.6, fontweight="bold")

axA.set_xlim(-0.7, 3.7)
axA.set_ylim(nrow - 0.5, -2.4)
axA.set_xticks([])
axA.set_yticks(range(nrow))
axA.set_yticklabels(sub.gene, fontsize=7.4, fontstyle="italic")
for s in ["top", "right", "bottom"]:
    axA.spines[s].set_visible(False)
axA.tick_params(axis="y", length=0, pad=2)
# module side bars
for ri, mod in enumerate(sub.module):
    axA.add_patch(Rectangle((-0.95, ri - .40), .16, .80, facecolor=COL[mod],
                            edgecolor="none", clip_on=False))
axA.text(-0.86, -1.4, "module", rotation=0, ha="center", va="bottom", fontsize=6.6, color="#4D4D4D")

# tier markers on the right
for ri, t in sub.iterrows():
    if str(t.tier).startswith("A_"):
        axA.text(3.62, ri, "\u25c0", va="center", ha="left", fontsize=5.4, color="#C0392B")
axA.text(3.62, -1.4, "Tier A", va="bottom", ha="left", fontsize=6.6, color="#C0392B")

leg = [Patch(facecolor=COL["SASP_inflammatory"], label="SASP / inflammatory"),
       Patch(facecolor=COL["ECM_fibrosis"], label="ECM / fibrosis"),
       Patch(facecolor="none", edgecolor="black", lw=1.0, label="padj < 0.05 (FAP)"),
       Patch(facecolor="#EFEFEF", edgecolor="white", label="not quantified")]
axA.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5, -0.015),
           ncol=2, frameon=False, fontsize=6.6, handlelength=1.2, columnspacing=1.2)

# ---------------------------------------------------------------- panel B
axB = fig.add_subplot(gs[0, 1])
fig.text(0.640, 0.985, "(b) FAP programs on the\n     population plasma age axis",
         fontsize=9.6, fontweight="bold", va="top")
gg = g.copy()
gg["label"] = gg.set.str.replace("_", " ")
gg = gg.sort_values("NES")
colr = ["#8E44AD" if "Clock" in s else "#2E86C1" if "Muscle" in s
        else "#E67E22" if "SASP" in s else "#27AE60" if "ECM" in s else "#8C8C8C" for s in gg.set]
bars = axB.barh(range(len(gg)), gg.NES, color=colr, alpha=.88, height=.62)
for i, (nes, p, n) in enumerate(zip(gg.NES, gg.p, gg.n_set)):
    star = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
    axB.text(9.4, i, f"{star} (n={n})", va="center", ha="right", fontsize=6.8,
             color="#404040" if star == "ns" else "black")
axB.axvline(0, color="black", lw=.7)
axB.set_yticks(range(len(gg)))
axB.set_yticklabels(gg.label, fontsize=7.2)
axB.set_xlabel("NES along the plasma protein age axis", fontsize=8)
axB.set_xlim(-6.8, 9.6)
for s in ["top", "right"]:
    axB.spines[s].set_visible(False)
axB.text(0, -0.145, "2,000 permutations; 1,848 plasma proteins with an FAP measurement.\n"
                    "Plasma age coefficients from Lehallier et al., Nat Med 2019 (n = 4,263).\n"
                    "In (A): red = higher with age / up in aged FAP, blue = lower; dark outline =\n"
                    "padj (FAP) or q (plasma) < 0.05, column 4 coloured by -log10(q).",
        transform=axB.transAxes, fontsize=6.4, color="#404040", va="top")

for ext, dpi in [("png", 600)]:
    fig.savefig(os.path.join(FIG, f"FigS9_plasma_translatability.{ext}"), dpi=dpi,
                bbox_inches="tight", facecolor="white")
print("saved FigS9_plasma_translatability.png")
