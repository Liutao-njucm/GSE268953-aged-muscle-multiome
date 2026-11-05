# -*- coding: utf-8 -*-
"""59_plasma_biomarker.py - map the FAP aging signature onto a population plasma proteome.

Plasma source : Lehallier et al., Nat Med 2019 (SomaScan, 4,263 donors, 2,925 assays;
                linear model per protein, age + sex + subcohort; BH q-values).
Signature     : our GSE268953 FAP pseudobulk DEG (snRNA) and GTEx muscle log2FC.

Outputs
  results/GSE268953/plasma/plasma_protein_age.csv          per-protein plasma age effect
  results/GSE268953/plasma/plasma_signature_stats.csv      enrichment tests
  results/GSE268953/plasma/plasma_biomarker_candidates.csv tiered candidate list
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
OUT = os.path.join(RES, "plasma")
os.makedirs(OUT, exist_ok=True)
XL = os.path.join(ROOT, "data", "plasma", "lehallier2019_SupplementaryTables.xlsx")

# ------------------------------------------------------------------ plasma
lm = pd.read_excel(XL, sheet_name="ST4 Linear modeling - Human", header=2)
lm = lm.dropna(subset=["ID"]).copy()
lm = lm.rename(columns={"ID": "plasma_id", "p.Age": "plasma_p", "q.Age": "plasma_q",
                        "Coefficient.Age": "plasma_coef_per_year"})
lm["plasma_coef_per_decade"] = lm["plasma_coef_per_year"] * 10

nom = pd.read_excel(XL, sheet_name="ST1 Nomenclature 2,925 proteins", header=2)
nom = nom.dropna(subset=["ID"]).copy()
nom = nom.rename(columns={"ID": "plasma_id", "Target": "plasma_target",
                          "TargetFullName": "plasma_target_full", "UniProt": "uniprot",
                          "EntrezGeneSymbol": "gene"})
nom["gene"] = nom["gene"].astype(str).str.strip()

pl = lm.merge(nom[["plasma_id", "gene", "plasma_target", "uniprot"]], on="plasma_id", how="inner")
pl = pl[pl["gene"].notna() & (pl["gene"] != "nan") & (pl["gene"] != "")]
# one row per gene: keep the assay with the smallest age q (most informative)
pl = pl.sort_values(["gene", "plasma_q"]).drop_duplicates("gene", keep="first")
pl["plasma_age_up"] = (pl["plasma_q"] < 0.05) & (pl["plasma_coef_per_decade"] > 0)
pl["plasma_age_down"] = (pl["plasma_q"] < 0.05) & (pl["plasma_coef_per_decade"] < 0)
pl = pl.reset_index(drop=True)
pl.to_csv(os.path.join(OUT, "plasma_protein_age.csv"), index=False)
print(f"plasma assays mapped to gene: {len(pl)}; age-associated (q<0.05): "
      f"{((pl.plasma_q < 0.05)).sum()} (up {pl.plasma_age_up.sum()}, down {pl.plasma_age_down.sum()})")

# ------------------------------------------------------------------ our signature
fap = pd.read_csv(os.path.join(RES, "deg", "DEG_FAP.csv"))
fap = fap.rename(columns={"log2FoldChange": "fap_lfc", "padj": "fap_padj"})
fap["gene"] = fap["gene"].astype(str)
fap["fap_sig_up"] = (fap.fap_padj < 0.05) & (fap.fap_lfc > 1)
fap["fap_sig_dn"] = (fap.fap_padj < 0.05) & (fap.fap_lfc < -1)

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

# GTEx muscle log2FC (for the third evidence layer)
gt = pd.read_csv(os.path.join(RES, "gtex", "gtex_muscle_validation.csv"))
gt = gt.rename(columns={"lfc_aged_vs_young": "gtex_lfc", "p_linear_model": "gtex_p"})
gt["gene"] = gt["gene"].astype(str)

# ------------------------------------------------------------------ merged table
m = fap.merge(pl, on="gene", how="inner").merge(gt[["gene", "gtex_lfc", "gtex_p"]], on="gene", how="left")

## ---- test 1: are FAP up-regulated DEGs enriched for plasma age-up proteins?
def fisher(a, b, c, d):
    orr, p = stats.fisher_exact([[a, b], [c, d]])
    return orr, p

sig_up_in_plasma = m[m.fap_sig_up]
n_up_plasma_ageup = int(sig_up_in_plasma.plasma_age_up.sum())
n_up_total = len(sig_up_in_plasma)
bg = m[~m.fap_sig_up]
n_bg_plasma_ageup = int(bg.plasma_age_up.sum())
n_bg_total = len(bg)
orr1, p1 = fisher(n_up_plasma_ageup, n_up_total - n_up_plasma_ageup,
                  n_bg_plasma_ageup, n_bg_total - n_bg_plasma_ageup)

sig_dn_in_plasma = m[m.fap_sig_dn]
orr2, p2 = fisher(int(sig_dn_in_plasma.plasma_age_down.sum()),
                  len(sig_dn_in_plasma) - int(sig_dn_in_plasma.plasma_age_down.sum()),
                  int(bg.plasma_age_down.sum()), len(bg) - int(bg.plasma_age_down.sum()))

## ---- test 2: global correlation between FAP log2FC and plasma age coefficient
rho, prho = stats.spearmanr(m.fap_lfc, m.plasma_coef_per_decade)
## restricted to plasma age-associated proteins only
sub = m[m.plasma_q < 0.05]
rho2, prho2 = stats.spearmanr(sub.fap_lfc, sub.plasma_coef_per_decade)

## ---- test 3: module-level shift of plasma age coefficient (one-sided MWU)
mod_rows = []
for k, g in m.groupby("module"):
    if k == "other" or len(g) < 3:
        continue
    u, pu = stats.mannwhitneyu(g.plasma_coef_per_decade, m[m.module == "other"].plasma_coef_per_decade,
                               alternative="two-sided")
    mod_rows.append(dict(module=k, n=len(g),
                         mean_plasma_coef=g.plasma_coef_per_decade.mean(),
                         mean_fap_lfc=g.fap_lfc.mean(),
                         frac_plasma_age_up=g.plasma_age_up.mean(),
                         frac_same_dir=(np.sign(g.fap_lfc) == np.sign(g.plasma_coef_per_decade)).mean(),
                         mwu_p=pu))
mods_df = pd.DataFrame(mod_rows)

stats_df = pd.DataFrame([
    dict(test="FAP up-DEGs (padj<0.05, LFC>1) vs background: fraction plasma age-up",
         n_case=n_up_total, k_case=n_up_plasma_ageup,
         n_ctrl=n_bg_total, k_ctrl=n_bg_plasma_ageup, est=orr1, p=p1, note="Fisher OR"),
    dict(test="FAP down-DEGs vs background: fraction plasma age-down",
         n_case=len(sig_dn_in_plasma), k_case=int(sig_dn_in_plasma.plasma_age_down.sum()),
         n_ctrl=len(bg), k_ctrl=int(bg.plasma_age_down.sum()), est=orr2, p=p2, note="Fisher OR"),
    dict(test="Spearman rho: FAP log2FC vs plasma age coef (all overlapping)",
         n_case=len(m), k_case=np.nan, n_ctrl=np.nan, k_ctrl=np.nan, est=rho, p=prho, note="rho"),
    dict(test="Spearman rho: restricted to plasma age-associated (q<0.05)",
         n_case=len(sub), k_case=np.nan, n_ctrl=np.nan, k_ctrl=np.nan, est=rho2, p=prho2, note="rho"),
])
stats_df = pd.concat([stats_df, mods_df.rename(columns={"module": "test"})], ignore_index=True)
stats_df.to_csv(os.path.join(OUT, "plasma_signature_stats.csv"), index=False)
print("\n=== enrichment tests ===")
print(stats_df.head(4).to_string(index=False))
print("\n=== module-level plasma age shift ===")
print(mods_df.to_string(index=False))

# ------------------------------------------------------------------ tiered candidates
# Tier 1: quantified in plasma panel AND age-associated (q<0.05) AND FAP log2FC > 0
# Tier 2: quantified in plasma, FAP log2FC > 0, but plasma age effect not significant
cand = m[(m.fap_lfc > 0)].copy()
cand["tier"] = np.where(cand.plasma_age_up, "Tier1_plasma_age_up",
                np.where(cand.plasma_q < 0.05, "Tier1x_plasma_age_down", "Tier2_plasma_flat"))
cand = cand.sort_values(["tier", "plasma_q"])
keep = ["gene", "module", "fap_lfc", "fap_padj", "fap_sig_up",
        "gtex_lfc", "gtex_p", "plasma_coef_per_decade", "plasma_q",
        "plasma_target", "uniprot", "tier"]
cand[keep].to_csv(os.path.join(OUT, "plasma_biomarker_candidates.csv"), index=False)
print(f"\ncandidates written: {len(cand)} (Tier1 {sum(cand.tier=='Tier1_plasma_age_up')}, "
      f"Tier2 {sum(cand.tier=='Tier2_plasma_flat')})")
print("\n=== Tier-1 candidates (plasma-quantified, age-up, FAP-up) ===")
t1 = cand[cand.tier == "Tier1_plasma_age_up"][keep]
print(t1.to_string(index=False))
