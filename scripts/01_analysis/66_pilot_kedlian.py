## 66_pilot_kedlian.py --- Kedlian Nat Aging 2024 stromal/endothelial compartments
## Independent replication test of the clock repressor arm (intercostal muscle, 3' v2/v3, whole cells)
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from scipy import stats

REP = ["NR1D1", "NR1D2", "PER1", "PER3", "DBP", "HLF", "TEF"]
ACT = ["ARNTL", "NPAS2", "CLOCK", "PER2", "CRY1"]
AGE_MID = {"15-20": 17.5, "25-30": 27.5, "35-40": 37.5, "50-55": 52.5, "55-60": 57.5,
           "60-65": 62.5, "70-75": 72.5}

def pseudobulk(path, keep_l1):
    a = ad.read_h5ad(path)
    m = a.obs.annotation_level1.isin(keep_l1)
    a = a[m].to_memory()
    X = sp.csr_matrix(a.X)
    sym = pd.Series(a.var["SYMBOL"].astype(str).values)
    codes, uniq = pd.factorize(sym)
    obs = a.obs.copy()
    obs["age_mid"] = obs["Age_group"].map(AGE_MID)
    donors = obs["donor_id"].astype(str).values
    rows, meta = [], []
    for d in pd.unique(donors):
        sel = donors == d
        cnt_full = np.asarray(X[sel].sum(axis=0)).ravel()
        cnt = np.zeros(len(uniq)); np.add.at(cnt, codes, cnt_full)
        cpm = np.log1p(1e4 * cnt / cnt.sum())
        rows.append(cpm)
        meta.append(dict(donor=d, age_group=obs["Age_group"][sel].iloc[0],
                         age=obs["age_mid"][sel].iloc[0],
                         bin=obs["Age_bin"][sel].iloc[0],
                         sex=obs["sex"][sel].iloc[0],
                         n_cells=int(sel.sum()),
                         median_umi=float(np.median(obs["n_counts"][sel]))))
    M = pd.DataFrame(rows, columns=list(uniq), index=[m["donor"] for m in meta])
    return M, pd.DataFrame(meta)

def analyze(M, meta, tag):
    res = []
    for g in REP + ACT:
        if g not in M.columns:
            res.append(dict(gene=g, missing=1)); continue
        m_old = (meta.bin == "old").values
        m_you = (meta.bin == "young").values
        y = M[g].values[m_old]
        yo = M[g].values[m_you]
        lfc = np.mean(y) - np.mean(yo)   # log1p-CPM10k difference
        u, p_two = stats.mannwhitneyu(y, yo, alternative="two-sided")
        rho, p_rho = stats.spearmanr(meta.age.values, M[g].values)
        res.append(dict(gene=g, arm="repressor" if g in REP else "other",
                        lfc_old_vs_young=round(lfc, 3),
                        pct_old_pos=np.mean(y > 0) * 100, pct_young_pos=np.mean(yo > 0) * 100,
                        mwu_p=round(p_two, 4), spearman_age_rho=round(rho, 3), spearman_p=round(p_rho, 4),
                        missing=0))
    d = pd.DataFrame([r for r in res if "arm" in r])
    d.to_csv(f"results/GSE167186_pilot/kedlian_{tag}_pseudobulk_test.csv", index=False)
    # module level (z-score genes across donors first)
    zz = M[[g for g in REP if g in M.columns]].apply(lambda c: (c - c.mean()) / c.std())
    meta[f"repressor_mod"] = zz.mean(axis=1).values
    zz2 = M[[g for g in ACT if g in M.columns]].apply(lambda c: (c - c.mean()) / c.std())
    meta[f"other_mod"] = zz2.mean(axis=1).values
    mo = (meta.bin == "old").values; my = (meta.bin == "young").values
    old = meta["repressor_mod"].values[mo]; yo = meta["repressor_mod"].values[my]
    u, p = stats.mannwhitneyu(old, yo, alternative="less")   # one-sided: old < young
    rho, pr = stats.spearmanr(meta.age.values, meta["repressor_mod"])
    ao = meta["other_mod"].values[mo]; ay = meta["other_mod"].values[my]
    ua, pa = stats.mannwhitneyu(ao, ay, alternative="two-sided")
    print(f"\n===== {tag} | donors old={sum(meta.bin=='old')} young={sum(meta.bin=='young')} "
          f"cells/donor median={int(meta.n_cells.median())} medianUMI={int(meta.median_umi.median())}")
    print(d.to_string(index=False))
    print(f"\nREPRESSOR module: one-sided old<young MWU p = {p:.4f}; Spearman vs age rho={rho:.2f} p={pr:.4f}")
    print(f"ACTIVATOR module: two-sided MWU p = {pa:.4f} (expect ns)")
    print(meta[["donor", "age", "bin", "sex", "n_cells", "repressor_mod", "other_mod"]].round(2).to_string(index=False))
    meta.round(4).to_csv(f"results/GSE167186_pilot/kedlian_{tag}_donor_module.csv", index=False)
    return d, meta

import os
os.makedirs("results/GSE167186_pilot", exist_ok=True)
M, meta = pseudobulk("data/Kedlian2024/kedlian_fibroblast_schwann.h5ad", ["FB"])
analyze(M, meta, "FAP_FB")

M2, meta2 = pseudobulk("data/Kedlian2024/kedlian_endothelium_smc.h5ad", ["CapEC", "VenEC", "ArtEC", "LymphEC"])
if len(meta2) >= 6:
    analyze(M2, meta2, "EC_SMC")
else:
    print("EC/SMC donors too few:", len(meta2))
