# -*- coding: utf-8 -*-
"""47_export_600dpi_compose.py
600-dpi re-export of the five main figures, exactly mirroring the layout of
scripts/27_assemble_main_figures.py (auto-trim + row layout) but at 2x scale:
  canvas 6000 px @ 600 dpi  (same physical 10-inch width as the 300-dpi originals)
Outputs TIFF (LZW) + PNG into submission/Figures_600dpi/.
Originals in submission/Figures/ are NOT touched.
"""
import os
from PIL import Image, ImageDraw, ImageFont

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
PAN = os.path.join(RES, "panels_600dpi")
OUT = os.path.join(RES, "submission", "Figures_600dpi")
os.makedirs(OUT, exist_ok=True)

CANVAS_W = 6000   # 2x the original 3000 px @ 300 dpi
GAP = 140
MARGIN = 80
LABEL_H = 180

try:
    FONT = ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 144)
except Exception:
    FONT = ImageFont.load_default()


def load_trim(path):
    im = Image.open(path)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(bg, im)
    im = im.convert("RGB")
    g = im.convert("L")
    bbox = g.point(lambda p: 0 if p > 252 else 255).getbbox()
    if bbox:
        pad = 60
        l, t, r, b = bbox
        l = max(0, l - pad); t = max(0, t - pad)
        r = min(im.width, r + pad); b = min(im.height, b + pad)
        im = im.crop((l, t, r, b))
    return im


def compose(rows, out_base):
    trimmed = [[(load_trim(p), 1) for p in row] for row in rows]
    nat = []
    for row in trimmed:
        ratios = [im.width / im.height for im, _ in row]
        h = CANVAS_W / sum(ratios)
        nat.append([(im, im.width / im.height * h, h) for im, _ in row])
    heights = [max(h for _, _, h in items) for items in nat]
    total_h = sum(heights) + GAP * (len(rows) - 1) + 2 * MARGIN + LABEL_H * len(rows)
    canvas = Image.new("RGB", (CANVAS_W, int(total_h)), "white")
    draw = ImageDraw.Draw(canvas)
    y = MARGIN
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    li = 0
    for items in nat:
        row_h = int(max(h for _, _, h in items))
        wsum = sum(w for _, w, _ in items)
        if wsum > CANVAS_W - 2 * MARGIN:
            f = (CANVAS_W - 2 * MARGIN - GAP * (len(items) - 1)) / wsum
        else:
            f = 1.0
        widths = [int(w * f) for _, w, _ in items]
        hs = [int(h * f) for _, _, h in items]
        row_w = sum(widths) + GAP * (len(items) - 1)
        x = (CANVAS_W - row_w) // 2
        for (im, _, _), w, h in zip(items, widths, hs):
            lab = letters[li]; li += 1
            draw.text((x + 8, y), lab, fill="black", font=FONT)
            img_y = y + LABEL_H
            r = im.resize((w, h), Image.LANCZOS)
            canvas.paste(r, (x, img_y))
            draw.rectangle([x - 1, img_y - 1, x + w, img_y + h], outline=(200, 200, 200))
            x += w + GAP
        y += LABEL_H + row_h + GAP
    tiff = os.path.join(OUT, out_base + ".tiff")
    canvas.save(tiff, dpi=(600, 600), compression="tiff_lzw")
    png = os.path.join(OUT, out_base + ".png")
    canvas.save(png, dpi=(600, 600))
    print(f"wrote {tiff}  {os.path.getsize(tiff)/1e6:.1f} MB | "
          f"{png}  {os.path.getsize(png)/1e6:.1f} MB  size={canvas.size}")


P = lambda name: os.path.join(PAN, name)

# ---------------- Fig 1: atlas ----------------
compose([
    [P("umap_celltype_600.png"), P("umap_group_600.png")],
    [P("celltype_proportion_stats_600.png")],
], "Figure1")

# ---------------- Fig 2: chromatin ----------------
compose([
    [P("p_DA_counts_600.png"), P("aging_TF_heatmap_600.png")],
    [P("chromatin_aging_score_by_celltype_600.png")],
], "Figure2")

# ---------------- Fig 3: transcriptome ----------------
compose([
    [P("p_DEG_counts_600.png")],
    [P("pathway_heatmap_600.png")],
    [P("FAP_clock_per_donor_600.png")],
], "Figure3")

# ---------------- Fig 4: coupling ----------------
compose([
    [P("coupling_FAP_600.png"), P("coupling_class_fraction_600.png")],
    [P("clock_TF_decoupling_FAP_600.png"), P("p_regchain_600.png")],
], "Figure4")

# ---------------- Fig 5: validation + model ----------------
compose([
    [P("gtex_concordance_by_module_600.png"), P("fig5_proteome_panel_600.png")],
    [P("clock_loci_tracks_600.png")],
    [P("p_model_600.png")],
], "Figure5")

print("ALL 600-DPI FIGURES DONE")
