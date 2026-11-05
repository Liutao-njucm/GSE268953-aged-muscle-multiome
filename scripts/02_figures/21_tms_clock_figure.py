# 21_tms_clock_figure.py — TMS mouse limb-muscle clock figure
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31/data/tms"
OUT  = r"D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"

r = pd.read_csv(f"{OUT}/tms_mouse_muscle_clock.csv", index_col=0)

# rebuild expression matrix for the plot (small: only needed genes x 54 samples)
import gzip, csv
meta = pd.read_csv(f"{ROOT}/metadata.csv.gz")
mus = meta[meta["source name"].str.startswith("Limb_Muscle")].copy()
mus["age"] = mus["characteristics: age"].astype(int)
mus["sex"] = mus["characteristics: sex"]
want = ["Nr1d1","Nr1d2","Per1","Per2","Per3","Dbp","Hlf","Tef","Arntl","Npas2",
        "Clock","Cry1","Cry2"]
cols = {}
with gzip.open(f"{ROOT}/counts.csv.gz","rt") as fh:
    header = fh.readline().rstrip("\n").split(",")
for c in header[1:]:
    sid = c.split(".gencode")[0]
    if sid in set(mus["Sample name"]): cols[c] = sid
idx = {g: i for i, g in enumerate(want)}
M = {g: {} for g in want}
with gzip.open(f"{ROOT}/counts.csv.gz","rt") as fh:
    rd = csv.reader(fh); next(rd)
    for rec in rd:
        if rec[0] not in idx: continue
        d = dict(zip(header[1:], rec[1:]))
        for c, sid in cols.items():
            M[rec[0]][sid] = int(d.get(c, 0))
df = pd.DataFrame(M)                      # rows = samples, columns = genes
df = df.reindex(columns=want).fillna(0)
lib = df.sum(axis=1)
log = np.log2(df.div(lib, axis=0) * 1e6 + 1).T   # rows = genes, cols = samples
m = mus.set_index("Sample name").loc[log.columns]
order = m.sort_values("age").index

fig, axes = plt.subplots(2, 7, figsize=(15.6, 5.4))
axes = axes.ravel()
## 断言：面板位必须容纳全部基因，否则 zip 会静默丢掉尾部基因
assert len(want) <= axes.size, f"{len(want)} genes but only {axes.size} panels"
for i, g in enumerate(want):
    ax = axes[i]
    y = log.loc[g, order]; a = m.loc[order, "age"]
    col = "#C0392B" if r.loc[g, "beta_age_per_month"] < 0 and r.loc[g, "p_age"] < 0.05 else \
          "#2E86C1" if r.loc[g, "p_age"] < 0.05 else "#7f7f7f"
    for s, c in [("m", "#1f77b4"), ("f", "#e377c2")]:
        sel = m.loc[order, "sex"] == s
        ax.scatter(a[sel], y[sel], s=18, c=c, alpha=0.75, edgecolors="none")
    # regression line over all
    k, b = np.polyfit(a, y, 1)
    xs = np.linspace(a.min(), a.max(), 50)
    ax.plot(xs, k*xs+b, color=col, lw=1.6)
    star = ""
    p = r.loc[g, "p_age"]
    if p < 0.001: star = "***"
    elif p < 0.01: star = "**"
    elif p < 0.05: star = "*"
    ax.set_title(f"{g} {star}", fontsize=10, fontweight="bold" if star else "normal")
    ax.tick_params(labelsize=8)
for ax in axes[len(want) + 1:]:      # 14 个面板位、13 个基因时仅剩的多余空位
    ax.axis("off")
# 图内图例：利用最后一个空面板位，读者无需翻图注即可读色
from matplotlib.lines import Line2D
leg_ax = axes[len(want)]
leg_ax.axis("off")
_handles = [
    Line2D([], [], marker="o", ls="none", color="#1f77b4", markersize=6, label="male"),
    Line2D([], [], marker="o", ls="none", color="#e377c2", markersize=6, label="female"),
    Line2D([], [], color="#2E86C1", lw=2, label="significant increase (p < 0.05)"),
    Line2D([], [], color="#C0392B", lw=2, label="significant decrease (p < 0.05)"),
    Line2D([], [], color="#7f7f7f", lw=2, label="not significant"),
    Line2D([], [], ls="none", label="age coefficient: * p<0.05  ** p<0.01  *** p<0.001"),
]
leg_ax.legend(handles=_handles, loc="center left", fontsize=9,
              frameon=False, title="Legend", title_fontsize=10)
for ax in axes[0::7]: ax.set_ylabel("log2 CPM", fontsize=9)     # 第一列
for ax in axes[7:]: ax.set_xlabel("age (months)", fontsize=9)   # 第二行
fig.suptitle("Tabula Muris Senis bulk limb muscle (n=54, ages 1-27 mo) — circadian clock across lifespan",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(f"{OUT}/figures/tms_mouse_clock_lifespan.png", dpi=220, bbox_inches="tight")
print("figure saved")
