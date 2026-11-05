# -*- coding: utf-8 -*-
"""47a_fig5_proteome_panel_600dpi.py
Re-render the Fig5 proteome panel (scripts/26) at 600 dpi into panels_600dpi/.
Does NOT overwrite figures/fig5_proteome_panel.png."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
OUT = os.path.join(RES, "panels_600dpi")
os.makedirs(OUT, exist_ok=True)

pan = pd.read_csv(os.path.join(RES, "proteome", "proteome_panel_genes.csv"))
sc = pd.read_csv(os.path.join(RES, "proteome", "proteome_module_scores_by_donor.csv"))

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.15, 1]})

# ---- left: ECM proteins ----
ax = axes[0]
de = pan[pan["panel"] == "sasp_ecm"].sort_values("log2FC_old_vs_young")
cols = ["#C0392B" if v > 0 else "#2E86C1" for v in de["log2FC_old_vs_young"]]
ax.barh(de["gene"], de["log2FC_old_vs_young"], color=cols, alpha=.88)
for y, (v, q) in enumerate(zip(de["log2FC_old_vs_young"], de["padj"])):
    star = "***" if q < 0.001 else ("**" if q < 0.01 else ("*" if q < 0.05 else ""))
    if star:
        ax.text(v + (0.03 if v > 0 else -0.03), y, star, va="center",
                ha="left" if v > 0 else "right", fontsize=9)
ax.axvline(0, c="k", lw=.8)
ax.set_xlabel("protein log2FC (old 65+ vs young 20-34)")
ax.set_title("ECM / matrisome proteins (n=15 detected)", fontsize=10, loc="left")
ax.tick_params(labelsize=8)

# ---- right: ECM score across age ----
ax = axes[1]
order = ["20-34", "35-49", "50-64", "65-79", "80+"]
data = [sc.loc[sc["age_group"] == g, "ecm_score"].dropna().values for g in order]
bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=.6, showfliers=False)
for patch in bp["boxes"]:
    patch.set_facecolor("#D5DBDB"); patch.set_edgecolor("#555")
rng = np.random.default_rng(7)
for i, dv in enumerate(data, start=1):
    ax.scatter(rng.normal(i, .07, len(dv)), dv, s=17, c="#7D3C98", alpha=.85, lw=0)
er, ep = stats.spearmanr(sc["age"], sc["ecm_score"])
try:
    import statsmodels.formula.api as smf
    m = smf.ols("ecm_score ~ age + C(sex)", data=sc.dropna(subset=["ecm_score"])).fit()
    ptxt = f"sex-adj. P = {m.pvalues['age']:.1e}"
except Exception:
    ptxt = ""
ax.axhline(0, c="k", lw=.5, ls=":")
ax.set_xlabel("donor age group"); ax.set_ylabel("ECM protein score (z)")
ax.set_title("ECM module rises with age (58 donors)", fontsize=10, loc="left")
ax.text(.03, .96, f"rho(age) = {er:+.2f}, P = {ep:.1e}\n{ptxt}",
        transform=ax.transAxes, va="top", fontsize=8.5)
ax.tick_params(labelsize=8)

fig.tight_layout()
out = os.path.join(OUT, "fig5_proteome_panel_600.png")
fig.savefig(out, dpi=600, bbox_inches="tight", facecolor="white")
print("wrote", out)
