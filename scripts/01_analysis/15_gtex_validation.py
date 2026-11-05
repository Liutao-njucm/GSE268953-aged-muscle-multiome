# -*- coding: utf-8 -*-
"""
Phase 2C-1: external validation in GTEx skeletal muscle (bulk RNA-seq, v10)
 - pull muscle samples (age/sex) and per-sample TPM for a clock/SASP/muscle panel
 - Young (20-39) vs Aged (60-79), linear model log2(TPM+1) ~ group + sex
 - cross-dataset concordance with our snRNA FAP log2FC
"""
import json, os, time, urllib.request, urllib.parse, csv, statistics

proj = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
out_dir = os.path.join(proj, "results/GSE268953/gtex")
os.makedirs(out_dir, exist_ok=True)
BASE = "https://gtexportal.org/api/v2"
DS = "gtex_v10"
TISSUE = "Muscle_Skeletal"

def get(path, **params):
    url = BASE + path + "?" + urllib.parse.urlencode(params)
    for a in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print("  retry", a, path, e, flush=True)
            time.sleep(3)
    return {}

# ---------- 1. muscle samples ----------
samples, page = [], 0
while True:
    d = get("/dataset/sample", datasetId=DS, tissueSiteDetailId=TISSUE,
            itemsPerPage=250, page=page)
    data = d.get("data", [])
    if not data:
        break
    samples += [(x["sampleId"], x.get("sex"), x.get("ageBracket"), x.get("dataType"))
                for x in data]
    page += 1
    if page >= d.get("paging_info", {}).get("numberOfPages", 0):
        break
print("muscle samples:", len(samples), flush=True)
rna = [s for s in samples if s[3] == "RNASEQ"]
rna = list({s[0]: s for s in rna}.values())
print("RNASEQ samples:", len(rna), flush=True)

# ---------- 2. gene panel ----------
genes = ["NR1D1","NR1D2","PER1","PER2","PER3","CRY1","CRY2","ARNTL","NPAS2","DBP","TEF",
         "HLF","NFIL3","CIART","BHLHE40","RORA","NR2F1","NR2F2","BHLHE41","SIRT1",
         "IL6","CXCL2","CCL2","MMP3","SERPINE1","IGFBP7","TGFB1","CSF1","IL1B","CXCL8",
         "CDKN1A","CDKN2A","GLB1","LMNB1",
         "MUSK","CRYAB","PFKM","PGK1","LMOD2","ATP5PD","STAC3","COL19A1","SAA2",
         "PDGFRA","DCN","LUM","POSTN","COL1A1","COL3A1","ADAMTS16","SEMA3A","PRG4",
         "ITGB8","MMP2","TIMP1","WISP2","CFD","ABCA8",
         "XIST","RPS4Y1","GAPDH","ACTB"]

def ensembl_id(sym):
    try:
        with urllib.request.urlopen(
            f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{sym}?content-type=application/json",
            timeout=60) as r:
            return json.loads(r.read().decode()).get("id")
    except Exception as e:
        print("  ensembl fail", sym, e, flush=True)
        return None

def gtex_gencode(ensg):
    d = get("/reference/gene", geneId=ensg, gencodeVersion="v39")
    data = d.get("data", [])
    return data[0]["gencodeId"] if data else None

def muscle_expr(gencode_id):
    d = get("/expression/geneExpression", gencodeId=gencode_id, datasetId=DS,
            tissueSiteDetailId=TISSUE, itemsPerPage=1)
    data = d.get("data", [])
    if not data:
        return None
    for blk in data:
        if blk.get("tissueSiteDetailId") == TISSUE:
            return blk.get("data")
    return data[0].get("data")

expr, meta = {}, {}
for i, g in enumerate(genes):
    e = ensembl_id(g)
    if not e:
        continue
    gid = gtex_gencode(e)
    if not gid:
        print("  no gtex id", g, flush=True); continue
    v = muscle_expr(gid)
    if not v:
        print("  no expr", g, flush=True); continue
    expr[g] = v
    meta[g] = (e, gid)
    print(f"{i+1}/{len(genes)} {g} -> {gid} n={len(v)}", flush=True)
    time.sleep(0.2)

# ---------- 3. check array ordering with sex-linked genes ----------
n_rna = len(rna)
print("\narray length:", {k: len(v) for k, v in list(expr.items())[:3]}, "rna samples:", n_rna)

sex_order = [s[1] for s in rna]
def sex_check(g):
    v = expr.get(g)
    if not v or len(v) != len(sex_order):
        return None
    fem = [x for x, s in zip(v, sex_order) if s == "female"]
    mal = [x for x, s in zip(v, sex_order) if s == "male"]
    if not fem or not mal:
        return None
    return statistics.mean(fem), statistics.mean(mal)

print("XIST female/male mean:", sex_check("XIST"))
print("RPS4Y1 female/male mean:", sex_check("RPS4Y1"))

# ---------- 4. statistics ----------
AGE_MID = {"20-29": 25, "30-39": 35, "40-49": 45, "50-59": 55, "60-69": 65, "70-79": 75}

def linreg(y, xs):
    """OLS with intercept + one binary covariate; returns (b_group, p_group) via t-test."""
    n = len(y)
    if n < 6: return None, None
    X = [[1.0] + list(x) for x in xs]
    k = len(X[0])
    # normal equations
    XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
    # gaussian elimination
    import copy
    A = [row[:] + [Xty[r]] for r, row in enumerate(XtX)]
    for c in range(k):
        piv = max(range(c, k), key=lambda r: abs(A[r][c]))
        A[c], A[piv] = A[piv], A[c]
        if abs(A[c][c]) < 1e-12: return None, None
        for r in range(k):
            if r != c:
                f = A[r][c]/A[c][c]
                for cc in range(c, k+1):
                    A[r][cc] -= f*A[c][cc]
    beta = [A[i][k]/A[i][i] for i in range(k)]
    # residual variance
    res = []
    for i in range(n):
        pred = sum(X[i][j]*beta[j] for j in range(k))
        res.append(y[i]-pred)
    s2 = sum(r*r for r in res)/(n-k)
    # se for beta1
    try:
        # invert XtX (k x k) via cofactor for k=2
        a, b, c_, d = XtX[0][0], XtX[0][1], XtX[1][0], XtX[1][1]
        det = a*d - b*c_
        inv11 = a/det
    except Exception:
        return beta[1], None
    se = (s2*inv11)**0.5
    if se == 0: return beta[1], None
    t = beta[1]/se
    # two-sided p from t with df = n-k (normal approx fine for n>50)
    import math
    p = math.erfc(abs(t)/math.sqrt(2))
    return beta[1], p

import math
rows = []
for g, v in expr.items():
    if len(v) != len(rna):
        continue
    ys, xs = [], []
    for val, s in zip(v, rna):
        ag = s[2]
        if ag in ("20-29", "30-39"):
            grp = 0
        elif ag in ("60-69", "70-79"):
            grp = 1
        else:
            continue
        ys.append(math.log2(val + 1))
        xs.append([grp, 1.0 if s[1] == "male" else 0.0])
    if len(ys) < 20:
        continue
    b, p = linreg(ys, xs)
    yng = [y for y, x in zip(ys, xs) if x[0] == 0]
    agd = [y for y, x in zip(ys, xs) if x[0] == 1]
    rows.append(dict(gene=g, n_young=len(yng), n_aged=len(agd),
                     mean_young_log2tpm=round(statistics.mean(yng), 3),
                     mean_aged_log2tpm=round(statistics.mean(agd), 3),
                     lfc_aged_vs_young=round(b, 3) if b is not None else "",
                     p_linear_model=f"{p:.3g}" if p is not None else ""))

rows.sort(key=lambda r: (r["lfc_aged_vs_young"] == "", r["lfc_aged_vs_young"] if isinstance(r["lfc_aged_vs_young"], float) else 0))
with open(os.path.join(out_dir, "gtex_muscle_validation.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

with open(os.path.join(out_dir, "gtex_muscle_samples.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["sampleId", "sex", "ageBracket", "dataType"]); w.writerows(rna)

print("\n=== GTEx skeletal muscle, Young(20-39) vs Aged(60-79) ===")
for r in rows:
    print(f"{r['gene']:10s} n={r['n_young']}/{r['n_aged']:3d}  log2FC={r['lfc_aged_vs_young']:>7}  p={r['p_linear_model']}")
print("\nsaved ->", out_dir)
