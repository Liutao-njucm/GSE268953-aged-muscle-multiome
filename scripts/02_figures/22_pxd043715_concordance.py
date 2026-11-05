# 22_pxd043715_concordance.py — protein-level concordance (mouse gastrocnemius old vs young) with human snRNA modules
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31/data/tms"
OUT  = r"D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"

## --- mouse proteome (PXD043715) ---
v = list(pd.read_excel(f"{ROOT}/pride_4type.xlsx", sheet_name=None).values())[0]
v["GN"] = v["Description"].str.extract(r"GN=([^ ]+)")
m = v[["GN", "Abundance Ratio: (Old) / (Young)",
       "Abundances (Grouped): Young", "Abundances (Grouped): Old"]].dropna(subset=["GN"])
m = m.rename(columns={"GN": "gene"})
m["log2FC_protein_old_vs_young"] = np.log2(m["Abundance Ratio: (Old) / (Young)"].astype(float))
m = m[m["log2FC_protein_old_vs_young"].abs() < 10]   # drop sentinel 100/0.01 imputations
print("proteins kept:", len(m))

## --- human modules ---
fap = pd.read_csv(f"{OUT}/deg/DEG_FAP.csv")
myo = pd.read_csv(f"{OUT}/deg/DEG_Myonuclei.csv")
def col(df, name):
    for c in df.columns:
        if c.lower() in ("gene", "gene_id", "x", "row"): return c
    return df.columns[0]
fapg = col(fap, "gene"); myog = col(myo, "gene")
fap_lfc = dict(zip(fap[fapg].astype(str), fap["log2FoldChange"]))
myo_lfc = dict(zip(myo[myog].astype(str), myo["log2FoldChange"]))

## map mouse gene -> human ortholog by capitalisation
m["human"] = m["gene"].str.upper()
m["fap_lfc"] = m["human"].map(fap_lfc)
m["myo_lfc"] = m["human"].map(myo_lfc)
m["in_FAP_DEG"] = m["fap_lfc"].notna()
m["in_Myo_DEG"] = m["myo_lfc"].notna()
m.to_csv(f"{OUT}/pxd043715_protein_concordance.csv", index=False)

ov = m[m["fap_lfc"].notna() | m["myo_lfc"].notna()]
print("\nOverlapping with human DEGs:", len(ov))
print(ov[["gene","human","log2FC_protein_old_vs_young","fap_lfc","myo_lfc"]].round(2).to_string(index=False))

for lab, hcol in [("FAP", "fap_lfc"), ("Myonuclei", "myo_lfc")]:
    s = m[m[hcol].notna()]
    if len(s) >= 5:
        rho, p = stats.spearmanr(s[hcol], s["log2FC_protein_old_vs_young"])
        agree = (np.sign(s[hcol]) == np.sign(s["log2FC_protein_old_vs_young"])).mean()
        print(f"{lab}: n={len(s)} Spearman rho={rho:.2f} p={p:.3f} sign-agreement={agree:.0%}")

## --- figure ---
def smart_labels(ax, xs, ys, texts, fontsize=6.5):
    """碰撞感知标签：依次尝试候选偏移，选第一个与已放置标签无重叠的位置；
    靠近右缘的点优先把文字放到点左侧，避免跨入右子图。"""
    fig_ = ax.figure
    fig_.canvas.draw()
    rend = fig_.canvas.get_renderer()
    x0, x1 = ax.get_xlim()
    placed = []
    cands = [(4, 3, "left", "bottom"), (4, -3, "left", "top"),
             (-4, 3, "right", "bottom"), (-4, -3, "right", "top"),
             (0, 7, "center", "bottom"), (0, -8, "center", "top"),
             (10, 0, "left", "center"), (-10, 0, "right", "center"),
             (7, 9, "left", "bottom"), (-7, 9, "right", "bottom"),
             (7, -10, "left", "top"), (-7, -10, "right", "top")]
    for x, y, t in zip(xs, ys, texts):
        order = list(cands)
        if x > x0 + 0.75 * (x1 - x0):          # 右缘：先试左侧
            order.sort(key=lambda c: c[2] != "right")
        chosen = None
        for dx, dy, ha, va in order:
            a = ax.annotate(t, (x, y), fontsize=fontsize, xytext=(dx, dy),
                            textcoords="offset points", ha=ha, va=va, zorder=5)
            bb = a.get_window_extent(renderer=rend)
            if not any(bb.overlaps(pb) for pb in placed):
                chosen = a
                break
            a.remove()
        if chosen is None:                      # 全部候选都撞：向上逐级堆叠找空位
            for k in range(1, 9):
                a = ax.annotate(t, (x, y), fontsize=fontsize, xytext=(0, 7 + k * 9),
                                textcoords="offset points", ha="center", va="bottom", zorder=5)
                bb = a.get_window_extent(renderer=rend)
                if not any(bb.overlaps(pb) for pb in placed):
                    chosen = a
                    break
                a.remove()
        if chosen is None:                      # 兜底（理论不达）
            chosen = ax.annotate(t, (x, y), fontsize=fontsize, xytext=(4, 3),
                                 textcoords="offset points", zorder=5)
        placed.append(chosen.get_window_extent(renderer=rend))

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
ax = axes[0]
s = m[m["fap_lfc"].notna()]
ax.scatter(s["fap_lfc"], s["log2FC_protein_old_vs_young"], c="#C0392B", s=40, alpha=0.8)
smart_labels(ax, s["fap_lfc"].values, s["log2FC_protein_old_vs_young"].values, s["gene"].values)
ax.axhline(0, color="grey", lw=0.8, ls="--"); ax.axvline(0, color="grey", lw=0.8, ls="--")
ax.set_xlabel("human FAP RNA log2FC (Aged vs Young)")
ax.set_ylabel("mouse gastrocnemius protein log2FC (Old/Young)")
ax.set_title("FAP-module proteins (PXD043715)", fontsize=10, fontweight="bold")
ax = axes[1]
s = m[m["myo_lfc"].notna()]
ax.scatter(s["myo_lfc"], s["log2FC_protein_old_vs_young"], c="#2E86C1", s=40, alpha=0.8)
smart_labels(ax, s["myo_lfc"].values, s["log2FC_protein_old_vs_young"].values, s["gene"].values)
ax.axhline(0, color="grey", lw=0.8, ls="--"); ax.axvline(0, color="grey", lw=0.8, ls="--")
ax.set_xlabel("human Myonuclei RNA log2FC (Aged vs Young)")
ax.set_ylabel("mouse gastrocnemius protein log2FC (Old/Young)")
ax.set_title("Myonuclei-module proteins (PXD043715)", fontsize=10, fontweight="bold")
fig.suptitle("Cross-species protein-level concordance (PXD043715, n_old=16 vs n_young=16)", fontweight="bold", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(f"{OUT}/figures/pxd043715_protein_concordance.png", dpi=220, bbox_inches="tight")
print("figure saved")
