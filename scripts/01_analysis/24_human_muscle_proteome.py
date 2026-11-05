#!/usr/bin/env python
"""24: donor-level analysis of the Ubaida-Mohien et al. eLife 2019 human skeletal
muscle TMT proteome (58 donors, 20-34 -> 80+ yr), and concordance with our
GSE268953 snMultiome FAP signature.

Input : data/proteome/elife-49874-fig1-data3-v1.xlsx  (Table S2, 5,891 proteins x 66 channels)
        data/proteome/uniprot_entry_to_gene.csv       (entry name -> gene symbol)
Output: results/GSE268953/proteome/human_muscle_proteome_aging.csv
        results/GSE268953/proteome/proteome_vs_fap_concordance.csv
        results/GSE268953/figures/human_proteome_validation.png
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROT = os.path.join(BASE, "data", "proteome")
RES = os.path.join(BASE, "results", "GSE268953", "proteome")
FIG = os.path.join(BASE, "results", "GSE268953", "figures")
os.makedirs(RES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------- label sheet
lab = pd.read_excel(os.path.join(PROT, "elife-49874-fig1-data3-v1.xlsx"),
                    sheet_name="Sample-Labeling_Details", header=None)
hdr_row = 1                                     # 0 = title, 1 = column names
lab = lab.iloc[hdr_row + 1:, :6].copy()
lab.columns = ["donor", "tmt", "channel", "age_group", "age", "sex"]
lab["donor"] = lab["donor"].astype(str).str.strip()
lab["tmt"] = lab["tmt"].astype(str).str.strip()
lab["channel"] = lab["channel"].astype(str).str.strip().str.replace(".0", "", regex=False)
lab["age"] = pd.to_numeric(lab["age"], errors="coerce")
lab = lab.dropna(subset=["age"])
print("donors with age:", len(lab))
print(lab["age_group"].value_counts().to_string())

# ---------------------------------------------------------------- protein table
raw = pd.read_excel(os.path.join(PROT, "elife-49874-fig1-data3-v1.xlsx"),
                    sheet_name="All-Proteins-5891-Quantified", header=None)
tmt_row, ch_row, first = 2, 3, 4
colkeys = []
for j in range(2, raw.shape[1]):
    colkeys.append((str(raw.iloc[tmt_row, j]).strip(), str(raw.iloc[ch_row, j]).strip().replace(".0", "")))
body = raw.iloc[first:, :].copy()
body = body[body.iloc[:, 1].notna()]
body.columns = ["num", "protein"] + [f"{t}|{c}" for t, c in colkeys]
body["protein"] = body["protein"].astype(str).str.strip()
for c in body.columns[2:]:
    body[c] = pd.to_numeric(body[c], errors="coerce")
print("proteins:", len(body), " channels:", body.shape[1] - 2)

# ---------------------------------------------------------------- long -> wide by donor
long = body.melt(id_vars=["protein"], value_vars=[c for c in body.columns[2:]],
                 var_name="key", value_name="log2r")
long[["tmt", "channel"]] = long["key"].str.split("|", expand=True)
long = long.merge(lab[["donor", "tmt", "channel", "age", "age_group", "sex"]],
                  on=["tmt", "channel"], how="inner")
print("matched observations:", len(long), " donors:", long["donor"].nunique())

# median-centre each TMT channel (removes batch / loading offsets; 131 = reference)
long["log2r"] = long["log2r"] - long.groupby("key")["log2r"].transform("median")

mat = long.pivot_table(index="protein", columns="donor", values="log2r", aggfunc="mean")
meta = lab.drop_duplicates("donor").set_index("donor").loc[[c for c in mat.columns]]
print("matrix:", mat.shape)

# ---------------------------------------------------------------- group stats
old_groups = ["65-79", "80+"]
is_old = meta["age_group"].isin(old_groups).values
is_young = (meta["age_group"] == "20-34").values
print("old n =", is_old.sum(), " young n =", is_young.sum(),
      " mid(35-64) n =", ((~is_old) & (~is_young)).sum())

vals = mat.values
old_v, yng_v = vals[:, is_old], vals[:, is_young]
with np.errstate(invalid="ignore"):
    lfc = np.nanmean(old_v, 1) - np.nanmean(yng_v, 1)
    t, p = stats.ttest_ind(old_v, yng_v, axis=1, equal_var=False, nan_policy="omit")
n_old = (~np.isnan(old_v)).sum(1)
n_yng = (~np.isnan(yng_v)).sum(1)

# age trend (Spearman) across all 58 donors
rho, prho = [], []
ages = meta["age"].values
for i in range(vals.shape[0]):
    x, y = ages, vals[i]
    ok = ~np.isnan(y)
    if ok.sum() >= 20 and np.nanstd(y) > 0:
        r, pp = stats.spearmanr(x[ok], y[ok])
    else:
        r, pp = np.nan, np.nan
    rho.append(r); prho.append(pp)

def bh(pv):
    pv = np.asarray(pv, float)
    ok = ~np.isnan(pv)
    out = np.full_like(pv, np.nan)
    q = pv[ok]
    n = len(q)
    order = np.argsort(q)
    ranked = q[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(n); tmp[order] = np.minimum(ranked, 1)
    out[ok] = tmp
    return out

stat = pd.DataFrame({
    "protein": mat.index,
    "log2FC_old_vs_young": lfc,
    "mean_old": np.nanmean(old_v, 1),
    "mean_young": np.nanmean(yng_v, 1),
    "n_old": n_old, "n_young": n_yng,
    "p_welch": p,
})
stat["padj"] = bh(stat["p_welch"].values)
stat["rho_age"] = rho
stat["p_rho"] = prho
stat["padj_rho"] = bh(np.array(prho, float))

mapping = pd.read_csv(os.path.join(PROT, "uniprot_entry_to_gene.csv"))
mapping["gene"] = mapping["gene"].fillna("").astype(str).str.split(" ").str[0]
stat = stat.merge(mapping[["entry", "gene", "accession"]],
                  left_on="protein", right_on="entry", how="left")
stat = stat.sort_values("p_welch")
stat.to_csv(os.path.join(RES, "human_muscle_proteome_aging.csv"), index=False)
print("\nsignificant proteins (padj<0.05):", int((stat["padj"] < 0.05).sum()),
      " up:", int(((stat["padj"] < 0.05) & (stat["log2FC_old_vs_young"] > 0)).sum()),
      " down:", int(((stat["padj"] < 0.05) & (stat["log2FC_old_vs_young"] < 0)).sum()))
print(stat.head(12)[["protein", "gene", "log2FC_old_vs_young", "p_welch", "padj", "rho_age"]].to_string(index=False))

# ---------------------------------------------------------------- panels of interest
CLOCK = ["NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "HLF", "TEF", "ARNTL", "NPAS2",
         "CLOCK", "CRY1", "CRY2", "CIART", "NFIL3", "RORA", "RORB", "RORC", "BHLHE40", "BHLHE41"]
SASP = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "LUM", "DCN",
        "POSTN", "FN1", "TIMP1", "TIMP3", "IGFBP7", "IGFBP3", "TGFB1", "TGFB2", "CTGF", "CCN2",
        "CCL2", "CXCL12", "IL6", "IL8", "CXCL1", "MMP2", "MMP14", "THBS1", "COMP", "ASPN",
        "SPARC", "SERPINE1", "IGF1", "MSTN", "GDF8", "VEGFA", "PDGFRA", "ABCA8", "MUSK",
        "CDKN2A", "CDKN1A", "TP53", "LMNA", "SIRT1", "SIRT6", "FOXO3", "MTOR"]
MYO = ["MYH1", "MYH2", "MYH4", "MYH7", "ACTA1", "ACTN2", "CKM", "LDHA", "LDHB", "PFKM",
       "PGK1", "PFKP", "ENO3", "GAPDH", "TPI1", "MYL1", "MYLPF", "TNNT3", "TNNI2", "ATP2A1"]

def pick(name, genes):
    d = stat[stat["gene"].isin(genes)].copy()
    d["panel"] = name
    return d

pan = pd.concat([pick("clock", CLOCK), pick("sasp_ecm", SASP), pick("myofiber_metab", MYO)])
pan = pan.sort_values(["panel", "p_welch"])
pan.to_csv(os.path.join(RES, "proteome_panel_genes.csv"), index=False)
print("\n=== panel coverage ===")
for g in ["clock", "sasp_ecm", "myofiber_metab"]:
    d = pan[pan["panel"] == g]
    print(f"{g}: covered {len(d)}/{len({'clock':CLOCK,'sasp_ecm':SASP,'myofiber_metab':MYO}[g])}"
          f"  sig(padj<0.05) {int((d['padj']<0.05).sum())}")
print(pan[pan["panel"].isin(["clock", "sasp_ecm"])]
      [["panel", "gene", "log2FC_old_vs_young", "p_welch", "padj", "rho_age"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ---------------------------------------------------------------- vs our DEGs (all cell types)
protcol = stat[["gene", "log2FC_old_vs_young", "p_welch", "padj", "rho_age"]]
conc_rows = []
comp_fap = None
comps = {}
for ct in ["FAP", "Myonuclei", "MuSC", "EC", "SMC_Pericyte", "T_NK_cells", "Macrophage"]:
    f = os.path.join(BASE, "results", "GSE268953", "deg", f"DEG_{ct}.csv")
    if not os.path.exists(f):
        continue
    d = pd.read_csv(f)
    d = d[["gene", "log2FoldChange", "padj"]].rename(columns={"log2FoldChange": "rna_lfc"})
    d["gene"] = d["gene"].astype(str).str.upper()
    c = d.merge(protcol, on="gene", how="inner").dropna(subset=["rna_lfc", "log2FC_old_vs_young"])
    if len(c) < 50:
        print(f"  {ct}: only {len(c)} overlapping genes, skipped"); continue
    r, pr = stats.spearmanr(c["rna_lfc"], c["log2FC_old_vs_young"])
    same = np.sign(c["rna_lfc"]) == np.sign(c["log2FC_old_vs_young"])
    # restrict to significant DEGs for a cleaner direction test
    sig = c[c["padj_x"] < 0.05] if "padj_x" in c.columns else c.iloc[0:0]
    conc_rows.append(dict(celltype=ct, n=len(c), rho=r, p_rho=pr,
                          same_dir_pct=100 * same.mean(),
                          n_sig_deg=len(sig),
                          same_dir_sig_pct=(100 * (np.sign(sig["rna_lfc"]) ==
                                                  np.sign(sig["log2FC_old_vs_young"])).mean()) if len(sig) else np.nan))
    print(f"  {ct}: n={len(c)}  rho={r:+.3f} (P={pr:.1e})  same-dir={100*same.mean():.1f}%")
    c["concordant"] = same
    comps[ct] = c.copy()
    if ct == "FAP":
        comp_fap = c.copy()

conc = pd.DataFrame(conc_rows).sort_values("rho", ascending=False)
conc.to_csv(os.path.join(RES, "proteome_concordance_by_celltype.csv"), index=False)
print("\n", conc.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

comp = comp_fap
r = float(conc.loc[conc.celltype == "FAP", "rho"].iloc[0])
pr = float(conc.loc[conc.celltype == "FAP", "p_rho"].iloc[0])
same = comp["concordant"]
comp.to_csv(os.path.join(RES, "proteome_vs_fap_concordance.csv"), index=False)

# ---------------------------------------------------------------- module scores
CLOCK_REP = ["NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "TEF", "HLF", "CIART", "NFIL3", "BHLHE41"]
ECM_SCORE = ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "TIMP3", "POSTN",
             "COMP", "LMNA", "DCN", "FN1", "ASPN", "SPARC", "IGFBP7"]


def score(genes):
    ps = stat.loc[stat["gene"].isin(genes), "protein"].tolist()
    ps = [p for p in ps if p in mat.index]
    if not ps:
        return None, []
    sub = mat.loc[ps]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-9, axis=0)
    return z.mean(axis=0), ps


rep, rep_ps = score(CLOCK_REP)
ecm, ecm_ps = score(ECM_SCORE)
n_clock_detected = int(stat["gene"].isin(CLOCK_REP).sum())
print(f"\nclock repressor proteins quantified: {len(rep_ps)}/{len(CLOCK_REP)} -> {rep_ps}")
print(f"ECM/Matrisome proteins quantified:  {len(ecm_ps)}/{len(ECM_SCORE)}")

sc = pd.DataFrame({"donor": mat.columns,
                   "ecm_score": ecm.values if ecm is not None else np.nan,
                   "clock_rep_score": rep.values if rep is not None else np.nan})
sc = sc.merge(meta.reset_index()[["donor", "age", "age_group", "sex"]], on="donor", how="left")
sc.to_csv(os.path.join(RES, "proteome_module_scores_by_donor.csv"), index=False)

if ecm is not None:
    er, ep = stats.spearmanr(sc["age"], sc["ecm_score"])
    tt = stats.ttest_ind(sc.loc[sc["age_group"].isin(old_groups), "ecm_score"],
                         sc.loc[sc["age_group"] == "20-34", "ecm_score"], equal_var=False)
    print(f"ECM protein score: rho(age)={er:+.3f} P={ep:.2e}; old-young diff={tt.statistic:+.2f} P={tt.pvalue:.2e}")
    try:
        import statsmodels.formula.api as smf
        m = smf.ols("ecm_score ~ age + C(sex)", data=sc.dropna(subset=["ecm_score"])).fit()
        print(f"  OLS ecm_score ~ age + sex: beta_age={m.params['age']:+.4f} P={m.pvalues['age']:.2e} "
              f"(n={int(m.nobs)}, R2={m.rsquared:.3f})")
    except Exception as e:
        print("  OLS skipped:", str(e)[:80])
else:
    er, ep = np.nan, np.nan
if rep is not None:
    rr, pp = stats.spearmanr(sc["age"], sc["clock_rep_score"])
    print(f"clock repressor score ({len(rep_ps)} protein): rho(age)={rr:+.3f} P={pp:.2e}")
else:
    rr, pp = np.nan, np.nan

# ---------------------------------------------------------------- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(19.5, 9.4))
gs = fig.add_gridspec(2, 4, hspace=0.45, wspace=0.38)

# A: volcano-style proteome aging
ax = fig.add_subplot(gs[0, 0])
d = stat.dropna(subset=["log2FC_old_vs_young", "p_welch"])
sig = d["padj"] < 0.05
ax.scatter(d.loc[~sig, "log2FC_old_vs_young"], -np.log10(d.loc[~sig, "p_welch"]),
           s=4, c="#B0B0B0", alpha=.45, lw=0)
ax.scatter(d.loc[sig, "log2FC_old_vs_young"], -np.log10(d.loc[sig, "p_welch"]),
           s=6, c="#C0392B", alpha=.8, lw=0)
ax.axvline(0, c="k", lw=.6, ls="--"); ax.axhline(-np.log10(0.05), c="k", lw=.6, ls=":")
ax.set_xlabel("log2FC (old vs young), protein"); ax.set_ylabel("-log10 P")
ax.set_title(f"a  Proteome-wide aging ({len(d):,} proteins)", fontsize=10, loc="left")
ax.text(.02, .96, f"{int(sig.sum())} proteins FDR<0.05\nyoung n={int(is_young.sum())}, old n={int(is_old.sum())}",
        transform=ax.transAxes, va="top", fontsize=8, color="#444")

# B: strongest age-associated proteins
ax = fig.add_subplot(gs[0, 1])
top = stat.dropna(subset=["padj"]).head(15).iloc[::-1]
cols = ["#C0392B" if v > 0 else "#2E86C1" for v in top["log2FC_old_vs_young"]]
ax.barh([g if isinstance(g, str) else p.split("_")[0] for g, p in zip(top["gene"], top["protein"])],
        top["log2FC_old_vs_young"], color=cols, alpha=.85)
ax.axvline(0, c="k", lw=.7); ax.set_xlabel("protein log2FC (old vs young)")
ax.set_title("b  Top age-associated proteins", fontsize=10, loc="left"); ax.tick_params(labelsize=7)

# C: SASP / ECM proteins
ax = fig.add_subplot(gs[0, 2])
ds = pan[pan["panel"] == "sasp_ecm"].sort_values("log2FC_old_vs_young")
cols = ["#C0392B" if v > 0 else "#2E86C1" for v in ds["log2FC_old_vs_young"]]
ax.barh(ds["gene"], ds["log2FC_old_vs_young"], color=cols, alpha=.85)
for y, (v, q) in enumerate(zip(ds["log2FC_old_vs_young"], ds["padj"])):
    if q < 0.05:
        ax.text(v + (0.02 if v > 0 else -0.02), y, "*", va="center",
                ha="left" if v > 0 else "right", fontsize=10)
ax.axvline(0, c="k", lw=.7); ax.set_xlabel("protein log2FC (old vs young)")
ax.set_title("c  SASP / ECM proteins", fontsize=10, loc="left"); ax.tick_params(labelsize=7)

# D: myofiber / metabolism proteins
ax = fig.add_subplot(gs[0, 3])
dm = pan[pan["panel"] == "myofiber_metab"].sort_values("log2FC_old_vs_young")
cols = ["#C0392B" if v > 0 else "#2E86C1" for v in dm["log2FC_old_vs_young"]]
ax.barh(dm["gene"], dm["log2FC_old_vs_young"], color=cols, alpha=.85)
ax.axvline(0, c="k", lw=.7); ax.set_xlabel("protein log2FC (old vs young)")
ax.set_title("d  Myofiber / glycolytic proteins", fontsize=10, loc="left"); ax.tick_params(labelsize=7)

# E: concordance with FAP RNA
ax = fig.add_subplot(gs[1, 0])
ax.scatter(comp["rna_lfc"], comp["log2FC_old_vs_young"], s=10,
           c=np.where(comp["concordant"], "#7D3C98", "#95A5A6"), alpha=.75, lw=0)
lim = [comp["rna_lfc"].quantile(.005), comp["rna_lfc"].quantile(.995)]
ax.plot(lim, lim, c="k", lw=.6, ls="--")
ax.axhline(0, c="k", lw=.5); ax.axvline(0, c="k", lw=.5)
ax.set_xlabel("FAP snRNA log2FC (aged/young)"); ax.set_ylabel("proteome log2FC")
ax.set_title(f"e  FAP signature vs bulk proteome (n={len(comp)})", fontsize=10, loc="left")
ax.text(.03, .96, f"Spearman rho={r:.2f} (P={pr:.1e})\nsame direction {100*same.mean():.0f}%",
        transform=ax.transAxes, va="top", fontsize=8.5)

# F: concordance by cell type
ax = fig.add_subplot(gs[1, 1])
cc = conc.sort_values("rho")
cols = ["#C0392B" if v > 0 else "#2E86C1" for v in cc["rho"]]
ax.barh(cc["celltype"], cc["rho"], color=cols, alpha=.85)
for y, (v, pv) in enumerate(zip(cc["rho"], cc["p_rho"])):
    star = "*" if pv < 0.05 else "ns"
    ax.text(v + (0.004 if v > 0 else -0.004), y, f"{v:+.2f} {star}", va="center",
            ha="left" if v > 0 else "right", fontsize=7.5)
ax.axvline(0, c="k", lw=.7)
ax.set_xlabel("Spearman rho (snRNA log2FC vs protein log2FC)")
ax.set_xlim(cc["rho"].min() - 0.045, cc["rho"].max() + 0.075)
ax.set_title("f  Which nuclear signature does bulk proteome track?", fontsize=10, loc="left")
ax.tick_params(labelsize=8)

# G: ECM protein score across age groups
ax = fig.add_subplot(gs[1, 2])
order = ["20-34", "35-49", "50-64", "65-79", "80+"]
data = [sc.loc[sc["age_group"] == g, "ecm_score"].dropna().values for g in order]
bp = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=.6, showfliers=False)
for patch in bp["boxes"]:
    patch.set_facecolor("#D5DBDB"); patch.set_edgecolor("#555")
for i, dvals in enumerate(data, start=1):
    ax.scatter(np.random.normal(i, .06, len(dvals)), dvals, s=16, c="#C0392B", alpha=.75, lw=0)
ax.axhline(0, c="k", lw=.5, ls=":")
ax.set_xlabel("age group"); ax.set_ylabel("ECM protein score (z)")
ax.set_title(f"g  ECM module, protein level ({len(ecm_ps)} proteins)", fontsize=10, loc="left")
ax.text(.03, .96, f"rho(age)={er:+.2f}, P={ep:.1e}", transform=ax.transAxes, va="top", fontsize=8.5)
ax.tick_params(labelsize=8)

# H: myonuclei signature concordance (the one that does track bulk)
ax = fig.add_subplot(gs[1, 3])
cm = comps.get("Myonuclei")
if cm is not None:
    ax.scatter(cm["rna_lfc"], cm["log2FC_old_vs_young"], s=10,
               c=np.where(cm["concordant"], "#E67E22", "#95A5A6"), alpha=.75, lw=0)
    lim = [cm["rna_lfc"].quantile(.005), cm["rna_lfc"].quantile(.995)]
    ax.plot(lim, lim, c="k", lw=.6, ls="--")
    ax.axhline(0, c="k", lw=.5); ax.axvline(0, c="k", lw=.5)
    mr = float(conc.loc[conc.celltype == "Myonuclei", "rho"].iloc[0])
    mp = float(conc.loc[conc.celltype == "Myonuclei", "p_rho"].iloc[0])
    ax.set_xlabel("Myonuclei snRNA log2FC"); ax.set_ylabel("proteome log2FC")
    ax.set_title(f"h  Myonuclear signature vs proteome (n={len(cm)})", fontsize=10, loc="left")
    ax.text(.03, .96, f"Spearman rho={mr:+.2f} (P={mp:.1e})", transform=ax.transAxes,
            va="top", fontsize=8.5)

fig.text(.5, .005,
         f"Note: only {len(rep_ps)}/{len(CLOCK_REP)} core clock proteins were quantified in this bulk proteome "
         f"({', '.join(rep_ps) if rep_ps else 'none'}), so the clock arm cannot be tested at protein level here.",
         ha="center", fontsize=8.5, style="italic", color="#555")

fig.suptitle("Human skeletal-muscle proteome across the lifespan (Ubaida-Mohien et al., eLife 2019; 58 donors, 5,891 proteins)",
             fontsize=12, y=.985)
out = os.path.join(FIG, "human_proteome_validation.png")
fig.savefig(out, dpi=250, bbox_inches="tight")
print("\nwrote", out)
