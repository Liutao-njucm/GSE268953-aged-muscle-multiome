"""
29_tf_regulon_activity.py

TF regulon activity from donor-level pseudobulk RNA, using CollecTRI (OmniPath)
and the Univariate Linear Model (ULM; decoupleR 2.x).

Rationale: the ATAC side of this study quantifies *motif accessibility* (chromVAR /
motif enrichment). This script adds the complementary *transcriptional output* layer:
does the regulon of a TF move in the same direction as its motif accessibility?
For the clock repressor arm this is the key test of the trans-acting model.

Input
  results/GSE268953/regulon/pb_logcpm_<celltype>.csv   (gene x donor log2 CPM)
  data/regulon/collectri_omnipath_2026.tsv             (TF -> target, signed)
Output
  results/GSE268953/regulon/regulon_activity_by_donor.csv
  results/GSE268953/regulon/regulon_activity_differential.csv
  results/GSE268953/regulon/regulon_panel_key_tfs.csv
"""
import os
import numpy as np
import pandas as pd
import decoupler as dc
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = "D:/workbuddy/outputs/2026-10-25-20-23-31"
REG = os.path.join(ROOT, "results/GSE268953/regulon")
NET = os.path.join(ROOT, "data/regulon/collectri_omnipath_2026.tsv")

CELLTYPES = ["FAP", "Myonuclei", "MuSC", "EC", "T_NK_cells", "Macrophage", "SMC_Pericyte"]

# ---------------------------------------------------------------- network
net_raw = pd.read_csv(NET, sep="\t")
print("CollecTRI raw edges:", len(net_raw))

def assign_mor(r):
    """+1 activation, -1 repression, NaN if direction unknown."""
    stim = str(r.get("consensus_stimulation", "")).lower() == "true"
    inhib = str(r.get("consensus_inhibition", "")).lower() == "true"
    if stim and not inhib:
        return 1
    if inhib and not stim:
        return -1
    return np.nan

net_raw["mor"] = net_raw.apply(assign_mor, axis=1)
net = net_raw.dropna(subset=["mor"])[["source_genesymbol", "target_genesymbol", "mor"]].copy()
net.columns = ["source", "target", "mor"]
net["mor"] = net["mor"].astype(int)
## resolve duplicated TF-target pairs: keep only unambiguous directions
amb = net.groupby(["source", "target"])["mor"].nunique()
net = net.set_index(["source", "target"]).loc[amb[amb == 1].index].reset_index()
net = net.drop_duplicates(subset=["source", "target"], keep="first")
net = net[net["source"] != net["target"]]
print(f"signed edges: {len(net)}  |  TFs: {net['source'].nunique()}  |  targets: {net['target'].nunique()}")

# ---------------------------------------------------------------- per cell type
donor_rows, diff_rows = [], []
for ct in CELLTYPES:
    f = os.path.join(REG, f"pb_logcpm_{ct}.csv")
    if not os.path.exists(f):
        print("missing", f); continue
    df = pd.read_csv(f).set_index("gene")
    df = df[[c for c in df.columns if c in df.columns]]
    donors = list(df.columns)
    grp = pd.Series(["Aged" if d.startswith("Aged") else "Young" for d in donors], index=donors)

    # per-gene z-score across donors (so ULM compares relative coordination, not baseline level)
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=1).replace(0, np.nan)
    z = df.sub(mu, axis=0).div(sd, axis=0).fillna(0.0)

    # keep TFs with enough measured targets in this cell type
    present = set(z.index)
    net_ct = net[net["target"].isin(present)].copy()

    est, pval = dc.mt.ulm(data=z.T, net=net_ct, tmin=5, verbose=False)
    # est: donors x TFs ; pval same shape (empty=True -> NaN where not tested)
    ntar = net_ct.groupby("source")["target"].nunique()

    est = est.loc[[d for d in donors if d in est.index]]
    for tf in est.columns:
        s = est[tf]
        y = s[grp.loc[s.index] == "Young"].dropna()
        a = s[grp.loc[s.index] == "Aged"].dropna()
        if len(y) < 3 or len(a) < 3:
            continue
        t, p = stats.ttest_ind(a, y, equal_var=False)
        for d in s.index:
            donor_rows.append({"celltype": ct, "donor": d, "group": grp[d],
                               "TF": tf, "score": float(s[d])})
        diff_rows.append({"celltype": ct, "TF": tf, "n_targets": int(ntar.get(tf, 0)),
                          "n_young": len(y), "n_aged": len(a),
                          "mean_young": y.mean(), "mean_aged": a.mean(),
                          "diff": a.mean() - y.mean(), "welch_t": t, "p": p})
    print(f"{ct:14s} donors={len(donors):2d}  TFs tested={est.shape[1]:4d}  "
          f"Young={int((grp=='Young').sum())} Aged={int((grp=='Aged').sum())}")

by_donor = pd.DataFrame(donor_rows)
dif = pd.DataFrame(diff_rows)
dif["padj"] = np.nan
for ct, idx in dif.groupby("celltype").groups.items():
    dif.loc[idx, "padj"] = multipletests(dif.loc[idx, "p"], method="fdr_bh")[1]
dif = dif.sort_values(["celltype", "p"])

by_donor.to_csv(os.path.join(REG, "regulon_activity_by_donor.csv"), index=False)
dif.to_csv(os.path.join(REG, "regulon_activity_differential.csv"), index=False)

# ---------------------------------------------------------------- key TF panel
CLOCK = ["NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "HLF", "TEF",
         "ARNTL", "ARNTL2", "CLOCK", "NPAS2", "BHLHE40", "BHLHE41", "CIART", "NFIL3"]
AP1 = ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "BATF", "BATF2", "BATF3", "MAFK", "MAFF"]
OTH = ["CEBPA", "CEBPB", "CEBPD", "CEBPE", "NFIA", "NFIB", "TEAD1", "TEAD4", "MEF2C",
       "PPARG", "EBF1", "EBF2", "RUNX2", "SMAD3", "STAT3", "NFKB1", "RELA", "YAP1", "WWTR1"]
panel = CLOCK + AP1 + OTH
key = dif[dif["TF"].isin(panel)].copy()
key["module"] = np.where(key["TF"].isin(CLOCK), "clock",
                         np.where(key["TF"].isin(AP1), "AP-1", "other"))
key = key[["celltype", "module", "TF", "n_targets", "mean_young", "mean_aged",
           "diff", "welch_t", "p", "padj"]]
key.to_csv(os.path.join(REG, "regulon_panel_key_tfs.csv"), index=False)

print("\n=== FAP, key TFs (sorted) ===")
show = key[(key["celltype"] == "FAP")].sort_values("p")
print(show.head(25).to_string(index=False, float_format=lambda v: f"{v:.3g}"))

print("\n=== clock TFs across cell types (padj<0.05) ===")
sig = key[(key["module"] == "clock") & (key["padj"] < 0.05)]
print(sig.sort_values(["celltype", "padj"]).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
print("\nDONE")
