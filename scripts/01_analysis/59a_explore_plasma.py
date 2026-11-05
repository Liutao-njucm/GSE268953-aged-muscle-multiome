# -*- coding: utf-8 -*-
"""59a_explore_plasma.py - explore intersection of the FAP SASP/ECM signature with the
Lehallier 2019 plasma proteome (SomaScan, 4,263 donors, 2,925 proteins)."""
import os
import numpy as np
import pandas as pd

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
XL = os.path.join(ROOT, "data", "plasma", "lehallier2019_SupplementaryTables.xlsx")

# ---- plasma: linear model results (age coefficient) ----
lm = pd.read_excel(XL, sheet_name="ST4 Linear modeling - Human", header=2)
lm = lm.dropna(subset=["ID"]).copy()
lm = lm[["ID", "p.Age", "q.Age", "Coefficient.Age"]].rename(
    columns={"ID": "plasma_id", "p.Age": "plasma_p", "q.Age": "plasma_q",
             "Coefficient.Age": "plasma_coef_per_year"})
lm["plasma_coef_per_decade"] = lm["plasma_coef_per_year"] * 10

nom = pd.read_excel(XL, sheet_name="ST1 Nomenclature 2,925 proteins", header=2)
nom = nom.dropna(subset=["ID"])[["ID", "Target", "EntrezGeneSymbol"]].copy()
nom = nom.rename(columns={"ID": "plasma_id", "Target": "plasma_target",
                          "EntrezGeneSymbol": "gene"})
nom["gene"] = nom["gene"].astype(str).str.strip()

pl = lm.merge(nom, on="plasma_id", how="left")
pl = pl[pl["gene"].notna() & (pl["gene"] != "nan")]
# one row per gene: keep the assay with the strongest age signal
pl["abs_q"] = pl["plasma_q"]
pl = pl.sort_values(["gene", "plasma_q"]).drop_duplicates("gene")
print("plasma proteins with gene symbol:", pl.shape[0])
print("age-associated in plasma (q<0.05):", (pl.plasma_q < 0.05).sum())

# ---- our FAP DEG ----
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
mod_of = {g: k for k, v in mods.items() for g in v}
fap["module"] = fap["gene"].map(mod_of).fillna("other")

m = fap.merge(pl, on="gene", how="inner")
print("\noverlap FAP DEGs x plasma proteins:", m.shape[0])
print(m["module"].value_counts().to_string())

print("\n--- module detail (mean FAP lfc / mean plasma coef per decade / same-dir %) ---")
for k, g in m.groupby("module"):
    same = (np.sign(g.fap_lfc) == np.sign(g.plasma_coef_per_decade)).mean() * 100
    print(f"{k:22s} n={len(g):4d}  fap={g.fap_lfc.mean():+.2f}  plasma={g.plasma_coef_per_decade.mean():+.3f}  same={same:.0f}%  "
          f"plasma_q<0.05 n={(g.plasma_q<0.05).sum()}")

print("\n--- SASP/ECM module genes detected in plasma ---")
se = m[m["module"].isin(["SASP_inflammatory", "ECM_fibrosis"])].sort_values("fap_lfc", ascending=False)
print(se[["gene","module","fap_lfc","fap_padj","plasma_coef_per_decade","plasma_q"]].to_string(index=False))
