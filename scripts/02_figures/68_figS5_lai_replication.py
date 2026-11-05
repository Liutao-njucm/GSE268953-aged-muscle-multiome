## 68_figS5_lai_replication.py --- independent replication of clock repressor arm in Lai Nat 2024
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

REP = ["NR1D1","NR1D2","PER1","PER3","DBP","HLF","TEF"]
ACT = ["ARNTL","NPAS2","CLOCK","PER2","CRY1"]
C_REP, C_ACT = "#B03A2E", "#2874A6"

M = pd.read_csv("results/GSE167186_pilot/lai_FAP_snRNA_donor_means.csv")
meta = pd.read_csv("results/GSE167186_pilot/lai_FAP_snRNA_donors.csv")
M.index = meta.donor.values
mod = M[REP].mean(axis=1).values
age = meta.age.values
sex = meta.sex.values

G = pd.read_csv("results/GSE167186_pilot/lai_FAP_all_clock_genes.csv")

fig = plt.figure(figsize=(7.2, 3.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.32,
                      left=0.09, right=0.97, top=0.80, bottom=0.17)

## (A) donor-level module vs age
axA = fig.add_subplot(gs[0, 0])
for s, c, mk in [("Male", "#1B2631", "o"), ("Female", "#E67E22", "^")]:
    m = sex == s
    axA.scatter(age[m], mod[m], s=42, c=c, marker=mk, label=f"{s} (n={m.sum()})",
                edgecolor="white", linewidth=0.6, zorder=3)
b, a0 = np.polyfit(age, mod, 1)
xs = np.linspace(age.min() - 2, age.max() + 2, 50)
axA.plot(xs, a0 + b * xs, color="#7B241C", lw=1.4, ls="--", zorder=2)
rho, pr = stats.spearmanr(age, mod)
mo = age >= 60; my = age < 50
u, pone = stats.mannwhitneyu(mod[mo], mod[my], alternative="less")
axA.set_title("(a) FAP clock repressor module vs age\n(Lai et al. Nature 2024, FAP snRNA, 13 donors)",
              loc="left", fontsize=8.6, fontweight="bold", pad=6)
axA.set_xlabel("Donor age (years)"); axA.set_ylabel("Repressor-arm module score\n(mean of scaled expression, 7 genes)")
axA.text(0.03, 0.05, f"Spearman rho = {rho:.2f} (p = {pr:.3f})\none-sided old<young p = {pone:.4f}",
         transform=axA.transAxes, fontsize=7.2, color="#7B241C")
axA.legend(fontsize=6.6, frameon=False, loc="upper right")
axA.tick_params(labelsize=7)

## (B) per-gene delta
axB = fig.add_subplot(gs[0, 1])
g = G.copy()
g = g.sort_values("arm", ascending=True)
ypos = np.arange(len(g))
colors = [C_REP if a == "repressor" else C_ACT for a in g.arm]
axB.barh(ypos, g.delta, color=colors, edgecolor="white", height=0.68)
axB.set_yticks(ypos); axB.set_yticklabels(g.gene, fontsize=7.4)
for i, (dl, p) in enumerate(zip(g.delta, g.mwu_p)):
    star = "**" if p < 0.01 else ("*" if p < 0.05 else "ns")
    off = 0.03 if dl >= 0 else -0.03
    axB.text(dl + off, i, star, va="center", ha="left" if dl >= 0 else "right",
             fontsize=7, color="#404040" if star == "ns" else "black")
axB.axvline(0, color="black", lw=0.8)
axB.set_title("(b) Old minus young, per gene\n(mean scaled expression; MWU: * p<0.05, ** p<0.01)",
              loc="left", fontsize=8.6, fontweight="bold", pad=6)
axB.set_xlabel("Delta (old - young)")
axB.tick_params(labelsize=7)
import matplotlib.patches as mpatches
axB.legend(handles=[mpatches.Patch(color=C_REP, label="repressor arm (module)"),
                    mpatches.Patch(color=C_ACT, label="other clock genes")],
           fontsize=6.6, frameon=False, loc="lower left", bbox_to_anchor=(0.02,0.02))

fig.suptitle("Independent replication of the FAP clock repressor-arm collapse\nin a public limb-muscle single-cell atlas",
             x=0.02, ha="left", fontsize=9.8, fontweight="bold", y=0.985)
fig.savefig("results/GSE268953/figures/FigS5_Lai_replication.png", dpi=300,
            facecolor="white", bbox_inches="tight")
print("saved FigS5")
