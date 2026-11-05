## 67_pilot_lai.py --- Lai et al. Nature 2024 (OMIX004308-02 stromal h5ad, limb muscle 15-99y)
## Usage: python scripts/67_pilot_lai.py inspect | run
import sys, numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
from scipy import stats

PATH = "data/Lai2024/OMIX004308-02_stromal.h5ad"
REP = ["NR1D1", "NR1D2", "PER1", "PER3", "DBP", "HLF", "TEF"]
ACT = ["ARNTL", "NPAS2", "CLOCK", "PER2", "CRY1"]
MOD = pd.read_csv("results/GSE268953/plasma/plasma_translatability.csv")
MODULES = {k: list(MOD[MOD.module == k].gene) for k in MOD.module.unique()}

def inspect():
    a = ad.read_h5ad(PATH, backed="r")
    print("shape", a.shape)
    print("obs cols:", list(a.obs.columns))
    for c in a.obs.columns:
        s = a.obs[c]
        if s.dtype.name in ("category", "object", "bool"):
            u = s.value_counts()
            print("---", c, f"({len(u)} levels)")
            print(u.head(25).to_string())
    a.file.close()

def run(l1_col, keep_vals, susp_col=None, susp_vals=None, tag="lai"):
    a = ad.read_h5ad(PATH, backed="r")
    m = a.obs[l1_col].isin(keep_vals)
    if susp_col:
        m &= a.obs[susp_col].isin(susp_vals)
    idx = np.where(m.values)[0]
    print(f"[{tag}] cells selected:", len(idx))
    sub = a[idx]
    X = sp.csr_matrix(sub.X)
    var = a.var
    symcol = "SYMBOL" if "SYMBOL" in var.columns else ("feature_name" if "feature_name" in var.columns else None)
    sym = pd.Series(var[symcol].astype(str).values) if symcol else pd.Series(var_names := [str(x) for x in a.var_names])
    codes, uniq = pd.factorize(sym)
    obs = sub.obs.copy()
    donors = obs["donor_id"].astype(str).values if "donor_id" in obs else obs.iloc[:, 0].astype(str).values
    rows, meta = [], []
    for d in pd.unique(donors):
        sel = donors == d
        cf = np.asarray(X[sel].sum(axis=0)).ravel()
        cnt = np.zeros(len(uniq)); np.add.at(cnt, codes, cf)
        rows.append(np.log1p(1e4 * cnt / cnt.sum()))
        meta.append(dict(donor=d, n_cells=int(sel.sum())))
    a.file.close()
    M = pd.DataFrame(rows, columns=list(uniq))
    meta = pd.DataFrame(meta)
    # try to derive age: find a column in original obs with numeric/age info per donor
    print("[{0}] donors: {1}".format(tag, len(meta)))
    return M, meta

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "inspect"
    if mode == "inspect":
        inspect()
    else:
        print("run mode: fill after inspecting obs")
