import gzip, json, time, csv, urllib.request, os

proj = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
feat = os.path.join(proj, "data/GSE268953/GSE268953_Young_1_features.tsv.gz")
out = os.path.join(proj, "data/annotation/ensembl_gene_coords.csv")

genes = []
with gzip.open(feat, "rt") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3 and p[2].strip() == "Gene Expression":
            genes.append(p[1].strip())
genes = sorted(set(genes))
print("genes in feature file:", len(genes))

def post(symbols):
    data = json.dumps({"symbols": symbols}).encode()
    req = urllib.request.Request(
        "https://rest.ensembl.org/lookup/symbol/homo_sapiens?content-type=application/json",
        data=data, headers={"Content-Type": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print("  retry", attempt, e)
            time.sleep(3)
    return {}

rows = {}
B = 900
for i in range(0, len(genes), B):
    chunk = genes[i:i+B]
    res = post(chunk)
    for sym, d in res.items():
        rows[sym] = (d.get("id"), d.get("seq_region_name"), d.get("start"), d.get("end"),
                     d.get("strand"), d.get("biotype"))
    print(f"{i+len(chunk)}/{len(genes)} done, mapped so far {len(rows)}", flush=True)
    time.sleep(0.4)

with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["symbol", "ensembl_id", "chr", "start", "end", "strand", "biotype"])
    for s in genes:
        if s in rows:
            w.writerow([s] + list(rows[s]))
print("mapped", len(rows), "/", len(genes), "->", out)
