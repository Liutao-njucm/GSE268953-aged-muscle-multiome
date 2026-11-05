#!/usr/bin/env python
"""24a: map UniProt entry names (e.g. 1433B_HUMAN) -> gene symbols.
Uses UniProt REST batch queries (id: query) because the stream endpoint times out here.
"""
import io, json, os, time, urllib.parse, urllib.request
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROT = os.path.join(BASE, "data", "proteome")
OUT = os.path.join(PROT, "uniprot_entry_to_gene.csv")


def fetch(ids):
    q = " OR ".join(f"id:{i}" for i in ids)
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q) +
           "&fields=accession,id,gene_primary&format=tsv&size=500")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                txt = r.read().decode("utf-8", "replace")
            rows = [ln.split("\t") for ln in txt.strip().split("\n")[1:] if ln.strip()]
            return [{"accession": a, "entry": e, "gene": g} for a, e, g in
                    (r + ["", "", ""])[:3] for r in [rows[0]]] if False else [
                {"accession": r[0], "entry": r[1], "gene": r[2] if len(r) > 2 else ""}
                for r in rows]
        except Exception as e:
            print("  retry", attempt, str(e)[:80]); time.sleep(3 * (attempt + 1))
    return []


def main():
    v = pd.read_excel(os.path.join(PROT, "elife-49874-fig1-data3-v1.xlsx"),
                      sheet_name="All-Proteins-5891-Quantified", header=None)
    entries = [str(x).strip() for x in v.iloc[4:, 1].tolist() if str(x).strip() not in ("", "nan")]
    entries = list(dict.fromkeys(entries))
    print("unique entry names:", len(entries))

    recs = []
    step = 80
    for i in range(0, len(entries), step):
        chunk = entries[i:i + step]
        recs += fetch(chunk)
        got = len({r["entry"] for r in recs})
        print(f"  chunk {i//step+1}/{(len(entries)+step-1)//step}  mapped={got}", flush=True)
        time.sleep(0.3)

    df = pd.DataFrame(recs).drop_duplicates(subset=["entry"])
    df.to_csv(OUT, index=False)
    print("mapped", len(df), "->", OUT)
    miss = [e for e in entries if e not in set(df["entry"])]
    print("unmapped", len(miss), miss[:10])


if __name__ == "__main__":
    main()
