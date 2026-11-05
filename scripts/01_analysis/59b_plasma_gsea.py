# -*- coding: utf-8 -*-
"""59b_plasma_gsea.py - gene-set enrichment of the FAP aging signature along the
plasma protein age axis (Lehallier 2019, 2,812 proteins ranked by age coefficient).

Rationale: 5-vs-5 snRNA has no power for single-gene calls, so we test whether the
FAP program as a whole sits high (or low) on the population plasma age axis.
"""
import os
import numpy as np
import pandas as pd

rng = np.random.default_rng(42)
ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
OUT = os.path.join(RES, "plasma")

pl = pd.read_csv(os.path.join(OUT, "plasma_protein_age.csv"))
fap = pd.read_csv(os.path.join(RES, "deg", "DEG_FAP.csv"))
fap = fap.rename(columns={"log2FoldChange": "fap_lfc", "padj": "fap_padj"})
fap["gene"] = fap["gene"].astype(str)

mods = {
    "SASP_inflammatory": ["IL6","CXCL2","CXCL8","CCL2","CSF1","TGFB1","SERPINE1","IGFBP7","MMP2",
                          "TIMP1","CDKN1A","CDKN2A","GLB1","LMNB1"],
    "ECM_fibrosis": ["COL1A1","COL3A1","DCN","LUM","POSTN","PDGFRA","PRG4","ADAMTS16","ITGB8",
                     "COL19A1","SAA2","ABCA8","MMP3","WISP2","CFD"],
    "Clock_repressive_arm": ["NR1D1","NR1D2","PER1","PER2","PER3","CRY1","CRY2","DBP","TEF","HLF",
                             "NR2F1","NR2F2","BHLHE40","BHLHE41","NFIL3","RORA","CIART"],
    "Muscle_fiber": ["PFKM","PGK1","LMOD2","STAC3","ATP5PD","CRYAB","MUSK"],
}
sig_up = set(fap[(fap.fap_padj < 0.05) & (fap.fap_lfc > 1)].gene)
sig_dn = set(fap[(fap.fap_padj < 0.05) & (fap.fap_lfc < -1)].gene)
top_up = set(fap.nlargest(150, "fap_lfc").gene)          # continuous alternative
top_dn = set(fap.nsmallest(150, "fap_lfc").gene)

# ranking universe: plasma proteins that are also detected in our FAP data
univ = pl[pl.gene.isin(set(fap.gene))].copy()
univ = univ.sort_values("plasma_coef_per_decade", ascending=False).reset_index(drop=True)
rank = univ.plasma_coef_per_decade.values
genes = univ.gene.tolist()
N = len(genes)
print(f"ranking universe: {N} plasma proteins with FAP measurement")


def es_score(hit_mask):
    """classic weighted Kolmogorov-Smirnov running sum (unweighted hits)."""
    nh = hit_mask.sum()
    if nh == 0 or nh == len(hit_mask):
        return 0.0
    steps = np.where(hit_mask, 1.0 / nh, -1.0 / (len(hit_mask) - nh))
    run = np.cumsum(steps)
    return run.max() if abs(run.max()) >= abs(run.min()) else run.min()


def gsea(set_genes, nperm=2000, label=""):
    s = set(set_genes) & set(genes)
    if len(s) < 4:
        return None
    mask = np.isin(genes, list(s))
    es = es_score(mask)
    null = np.empty(nperm)
    for i in range(nperm):
        pm = rng.permutation(mask)
        null[i] = es_score(pm)
    pos = null[null >= 0].std() if (null >= 0).any() else np.nan
    neg = null[null < 0].std() if (null < 0).any() else np.nan
    nes = es / pos if es >= 0 else es / abs(neg)
    p = (np.sum(null >= es) + 1) / (nperm + 1) if es >= 0 else (np.sum(null <= es) + 1) / (nperm + 1)
    # mean rank percentile of the set
    pct = np.mean([1 - i / (N - 1) for i in np.where(mask)[0]])
    return dict(set=label, n_set=len(s), ES=es, NES=nes, p=p, mean_rank_pct=pct,
                n_age_up_q05=int(univ.loc[mask, "plasma_q"].lt(0.05).sum()))


rows = []
for name, ss in [("FAP_up_DEG", sig_up), ("FAP_down_DEG", sig_dn),
                 ("FAP_top150_up", top_up), ("FAP_top150_down", top_dn),
                 ("SASP_inflammatory", mods["SASP_inflammatory"]),
                 ("ECM_fibrosis", mods["ECM_fibrosis"]),
                 ("Clock_repressive_arm", mods["Clock_repressive_arm"]),
                 ("Muscle_fiber", mods["Muscle_fiber"])]:
    r = gsea(ss, label=name)
    if r:
        rows.append(r)
        print(f"{name:22s} n={r['n_set']:3d}  ES={r['ES']:+.3f} NES={r['NES']:+.2f} "
              f"p={r['p']:.4f}  mean_rank_pct={r['mean_rank_pct']:.2f}")

g = pd.DataFrame(rows)
g.to_csv(os.path.join(OUT, "plasma_gsea_by_geneset.csv"), index=False)

# also: FAP lfc vs plasma age coef inside each module (direction agreement)
m = fap.merge(pl, on="gene", how="inner")
mod_of = {x: k for k, v in mods.items() for x in v}
m["module"] = m.gene.map(mod_of).fillna("other")
print("\nmodule direction agreement (FAP lfc vs plasma age coef):")
for k, d in m.groupby("module"):
    if k == "other":
        continue
    print(f"  {k:22s} n={len(d):2d} same_dir={np.mean(np.sign(d.fap_lfc)==np.sign(d.plasma_coef_per_decade)):.2f} "
          f"plasma_age_up={int(d.plasma_q.lt(0.05).mul(d.plasma_coef_per_decade.gt(0)).sum())}")
print("\nDONE")
