# 20_tms_clock_check.py — Tabula Muris Senis bulk RNA-seq, limb muscle across age
# Cross-species check of the clock repressor-arm collapse seen in aged human FAPs
import gzip, csv
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31/data/tms"
OUT  = r"D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"

meta = pd.read_csv(f"{ROOT}/metadata.csv.gz")
meta["sample_id"] = meta["Sample name"]
mus = meta[meta["source name"].str.startswith("Limb_Muscle")].copy()
mus["age"] = mus["characteristics: age"].astype(int)
mus["sex"] = mus["characteristics: sex"]
print("Limb_Muscle samples:", len(mus))
print(mus.groupby("age").size().to_string())

## load counts, keep muscle columns
cols_of_interest = {}
with gzip.open(f"{ROOT}/counts.csv.gz", "rt") as fh:
    header = fh.readline().rstrip("\n").split(",")
for c in header[1:]:
    sid = c.split(".gencode")[0]
    cols_of_interest[c] = sid
mus_ids = set(mus["sample_id"])
keep = [c for c, sid in cols_of_interest.items() if sid in mus_ids]
print("matching muscle count columns:", len(keep))

genes = []
rows = []
with gzip.open(f"{ROOT}/counts.csv.gz", "rt") as fh:
    rd = csv.reader(fh)
    next(rd)
    for rec in rd:
        vals = {c: int(v) for c, v in zip(header[1:], rec[1:]) if v}
        rows.append([rec[0]] + [vals.get(c, 0) for c in keep])

cnt = pd.DataFrame(rows).set_index(0)
cnt.columns = [c.split(".gencode")[0] for c in keep]
print("matrix:", cnt.shape)

## CPM + log2
lib = cnt.sum(axis=0)
cpm = cnt / lib * 1e6
log = np.log2(cpm + 1)
m = mus.set_index("sample_id").loc[cnt.columns]
assert (m.index == cnt.columns).all()

CLOCK = ["Nr1d1","Nr1d2","Per1","Per2","Per3","Dbp","Hlf","Tef","Arntl","Npas2",
         "Clock","Cry1","Cry2","Rora","Bhlhe41","E4bp4"]
CLOCK = list(dict.fromkeys(CLOCK))
CLOCK = [g for g in CLOCK if g in log.index]
SASP  = [g for g in ["Col1a1","Col1a2","Lum","Postn","Timp1","Ctgf","Il6","Ccl2",
                     "Vim","Pdgfra","Apoa1","Fbn1","Dcn"] if g in log.index]
MYO   = [g for g in ["Tnnt3","Tnnc2","Myh2","Actn2","Ckm","Tcap","Ldha","Pfkm"] if g in log.index]
print("clock genes found:", CLOCK)

res = []
for gene in CLOCK + SASP + MYO:
    y = log.loc[gene]
    # continuous age model (log2CPM ~ age + sex)
    X = np.column_stack([np.ones(len(m)), m["age"], (m["sex"] == "f").astype(float)])
    beta, *_ = np.linalg.lstsq(X, y.values, rcond=None)
    resid = y.values - X @ beta
    dof = len(m) - X.shape[1]
    s2 = resid @ resid / dof
    se_age = np.sqrt(s2 * np.linalg.pinv(X.T @ X)[1, 1])
    t_age = beta[1] / se_age
    p_age = 2 * stats.t.sf(abs(t_age), dof)
    # young (<=6m) vs old (>=24m) Welch
    yng = y[m["age"] <= 6]; old = y[m["age"] >= 24]
    tt = stats.ttest_ind(old, yng, equal_var=False)
    res.append(dict(gene=gene, module="Clock" if gene in CLOCK else ("SASP/ECM" if gene in SASP else "Myofibre"),
                    beta_age_per_month=beta[1], t_age=t_age, p_age=p_age,
                    young_mean=yng.mean(), old_mean=old.mean(),
                    log2FC_old_vs_young=old.mean() - yng.mean(), p_welch=tt.pvalue))
r = pd.DataFrame(res)
r["p_adj"] = pd.Series(stats.false_discovery_control(r["p_age"].values) if hasattr(stats, "false_discovery_control") else r["p_age"],
                       index=r.index)
try:
    from statsmodels.stats.multitest import multipletests
    r["p_adj_age"] = multipletests(r["p_age"], method="fdr_bh")[1]
    r["p_adj_welch"] = multipletests(r["p_welch"], method="fdr_bh")[1]
except ImportError:
    r["p_adj_age"] = r["p_adj"]; r["p_adj_welch"] = r["p_welch"]

r.to_csv(f"{OUT}/tms_mouse_muscle_clock.csv", index=False)
print(r.round(4).to_string(index=False))
