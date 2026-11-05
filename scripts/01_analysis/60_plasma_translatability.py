# -*- coding: utf-8 -*-
"""60_plasma_translatability.py - build the plasma translatability table for the
FAP SASP/ECM program.

For every gene in our SASP / ECM / clock / muscle modules:
  * FAP snRNA log2FC / padj (this study)
  * GTEx skeletal muscle log2FC (this study)
  * SomaScan plasma age coefficient + q (Lehallier 2019, 4,263 donors)
  * HPA secretome location and measured blood concentration (panel-independent)
Then classify each gene into translational tiers.

Inputs : results/GSE268953/plasma/plasma_protein_age.csv
         results/GSE268953/deg/DEG_FAP.csv, gtex/gtex_muscle_validation.csv
         APIs: mygene.info (symbol -> ENSG), proteinatlas.org (ENSG -> blood data)
Outputs: results/GSE268953/plasma/plasma_translatability.csv
         results/GSE268953/plasma/hpa_blood_raw.json
"""
import json
import os
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
OUT = os.path.join(RES, "plasma")
CACHE = os.path.join(OUT, "hpa_blood_raw.json")

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
genes = list(mod_of)

# ------------------------------------------------------------------ ENSG via mygene
def mygene_ensg(symbols):
    url = "https://mygene.info/v3/query"
    body = urllib.parse.urlencode({
        "q": ",".join(symbols), "scopes": "symbol", "species": "human",
        "fields": "ensembl.gene,symbol", "size": 1}).encode()
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(req, timeout=60) as r:
        hits = json.load(r)
    out = {}
    for h in hits:
        if "symbol" not in h:
            continue
        ens = h.get("ensembl")
        if isinstance(ens, list):
            ens = ens[0]
        if isinstance(ens, dict) and ens.get("gene"):
            out[h["symbol"]] = ens["gene"]
    return out


if os.path.exists(CACHE):
    hpa = json.load(open(CACHE))
    print(f"HPA cache loaded: {len(hpa)} genes")
else:
    hpa = {}

missing = [g for g in genes if g not in hpa]
if missing:
    try:
        e2g = mygene_ensg([g for g in missing if g not in hpa][:1000])
    except Exception as e:
        print("mygene batch failed:", e)
        e2g = {}
    print(f"ENSG resolved: {len(e2g)}/{len(missing)}")
    BL_KEY = "Blood concentration - Conc. blood MS [pg/L]"
    for g in missing:
        ensg = e2g.get(g)
        if not ensg:
            hpa[g] = {"ensg": None}
            continue
        url = f"https://www.proteinatlas.org/{ensg}.json"
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                d = json.load(r)
            hpa[g] = {"ensg": ensg,
                      "secretome_location": d.get("Secretome location"),
                      "secretome_function": d.get("Secretome function"),
                      "blood_conc_MS_pg_L": d.get(BL_KEY),
                      "blood_conc_IM_pg_L": d.get("Blood concentration - Conc. blood IM [pg/L]"),
                      "protein_class": d.get("Protein class")}
        except Exception as e:
            hpa[g] = {"ensg": ensg, "error": str(e)}
        time.sleep(0.25)
    json.dump(hpa, open(CACHE, "w"), indent=1)
print(f"HPA queried: {sum(1 for v in hpa.values() if v.get('ensg'))} genes")

# ------------------------------------------------------------------ merge evidence
pl = pd.read_csv(os.path.join(OUT, "plasma_protein_age.csv"))
fap = pd.read_csv(os.path.join(RES, "deg", "DEG_FAP.csv")).rename(
    columns={"log2FoldChange": "fap_lfc", "padj": "fap_padj"})
gt = pd.read_csv(os.path.join(RES, "gtex", "gtex_muscle_validation.csv")).rename(
    columns={"lfc_aged_vs_young": "gtex_lfc", "p_linear_model": "gtex_p"})

def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


rows = []
for g in genes:
    h = hpa.get(g, {})
    f = fap[fap.gene == g]
    t = gt[gt.gene == g]
    p = pl[pl.gene == g]
    conc_pg_L = _num(h.get("blood_conc_MS_pg_L"))
    if np.isnan(conc_pg_L):
        conc_pg_L = _num(h.get("blood_conc_IM_pg_L"))
    rows.append(dict(
        gene=g, module=mod_of[g],
        fap_lfc=float(f.fap_lfc.iloc[0]) if len(f) else np.nan,
        fap_padj=float(f.fap_padj.iloc[0]) if len(f) else np.nan,
        gtex_lfc=float(t.gtex_lfc.iloc[0]) if len(t) else np.nan,
        ssm_quantified=(len(p) > 0),
        plasma_coef_per_decade=float(p.plasma_coef_per_decade.iloc[0]) if len(p) else np.nan,
        plasma_q=float(p.plasma_q.iloc[0]) if len(p) else np.nan,
        hpa_secretome=h.get("secretome_location"),
        hpa_blood_conc_mg_L=conc_pg_L / 1e9 if not np.isnan(conc_pg_L) else np.nan,
    ))
d = pd.DataFrame(rows)
d["plasma_detectable"] = (d.ssm_quantified
                          | d.hpa_blood_conc_mg_L.notna()
                          | d.hpa_secretome.fillna("").str.contains("blood"))
d["plasma_age_up"] = (d.plasma_q < 0.05) & (d.plasma_coef_per_decade > 0)
d["fap_up"] = d.fap_lfc > 0
d["fap_sig_up"] = (d.fap_padj < 0.05) & (d.fap_lfc > 1)

def tier(r):
    if not r.plasma_detectable:
        return "C_not_measurable"
    if r.plasma_age_up and r.fap_up:
        return "A_detectable+age_up+FAP_up"
    if r.plasma_age_up:
        return "B_age_up_but_FAP_down"
    return "B2_detectable_no_age_signal"

d["tier"] = d.apply(tier, axis=1)
d = d.sort_values(["module", "tier", "plasma_q"])
d.to_csv(os.path.join(OUT, "plasma_translatability.csv"), index=False)

pd.set_option("display.width", 200)
print("\n=== translatability table (module genes) ===")
print(d[["gene", "module", "fap_lfc", "fap_padj", "gtex_lfc", "ssm_quantified",
         "plasma_coef_per_decade", "plasma_q", "hpa_blood_conc_mg_L", "tier"]]
      .to_string(index=False))
print("\ntier counts:\n", d.tier.value_counts().to_string())
