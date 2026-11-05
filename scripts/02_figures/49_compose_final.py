# -*- coding: utf-8 -*-
"""49_compose_final.py — assemble the five main figures from panels_600dpi/
at two scales:
  scale 1 (3000 px @ 300 dpi) -> submission/Figures/Figure1..5.png  (submission)
                                 + figures_main/Fig1..Fig5.png      (working copies)
  scale 2 (6000 px @ 600 dpi) -> submission/Figures_600dpi/Figure1..5.tiff (LZW) + .png
Layout mirrors scripts/27_assemble_main_figures.py."""
import os
from PIL import Image, ImageDraw, ImageFont

ROOT = r"D:/workbuddy/outputs/2026-10-25-20-23-31"
RES = os.path.join(ROOT, "results", "GSE268953")
PAN = os.path.join(RES, "panels_600dpi")

FIGS = {
    "Figure1": [[ "umap_celltype_600.png", "umap_group_600.png"],
                ["celltype_proportion_stats_600.png"]],
    "Figure2": [[ "p_DA_counts_600.png", "aging_TF_heatmap_600.png"],
                ["chromatin_aging_score_by_celltype_600.png"]],
    "Figure3": [[ "p_DEG_counts_600.png"],
                ["pathway_heatmap_600.png"],
                ["FAP_clock_per_donor_600.png"]],
    "Figure4": [[ "coupling_FAP_600.png", "coupling_class_fraction_600.png"],
                ["clock_TF_decoupling_FAP_600.png", "p_regchain_600.png"]],
    "Figure5": [[ "gtex_concordance_by_module_600.png", "fig5_proteome_panel_600.png"],
                ["clock_loci_tracks_600.png"],
                ["p_model_600.png"]],
}


def make_font(size):
    try:
        return ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", size)
    except Exception:
        return ImageFont.load_default()


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
        pad = int(30 * 2)
        l, t, r, b = bbox
        l = max(0, l - pad); t = max(0, t - pad)
        r = min(im.width, r + pad); b = min(im.height, b + pad)
        im = im.crop((l, t, r, b))
    return im


def compose(rows, base, canvas_w, gap, margin, label_h, font, out_pngs, out_tiffs):
    trimmed = [[load_trim(os.path.join(PAN, p)) for p in row] for row in rows]
    nat = []
    for row in trimmed:
        ratios = [im.width / im.height for im in row]
        h = canvas_w / sum(ratios)
        nat.append([(im, im.width / im.height * h, h) for im in row])
    heights = [max(h for _, _, h in items) for items in nat]
    total_h = sum(heights) + gap * (len(rows) - 1) + 2 * margin + label_h * len(rows)
    canvas = Image.new("RGB", (canvas_w, int(total_h)), "white")
    draw = ImageDraw.Draw(canvas)
    y = margin
    letters = "abcdefghijklmnopqrstuvwxyz"  # 与正文/图注的小写面板字母体系一致
    li = 0
    for items in nat:
        row_h = int(max(h for _, _, h in items))
        wsum = sum(w for _, w, _ in items)
        if wsum > canvas_w - 2 * margin:
            f = (canvas_w - 2 * margin - gap * (len(items) - 1)) / wsum
        else:
            f = 1.0
        widths = [int(w * f) for _, w, _ in items]
        hs = [int(h * f) for _, _, h in items]
        row_w = sum(widths) + gap * (len(items) - 1)
        x = (canvas_w - row_w) // 2
        for (im, _, _), w, h in zip(items, widths, hs):
            draw.text((x + 8, y), letters[li], fill="black", font=font); li += 1
            img_y = y + label_h
            canvas.paste(im.resize((w, h), Image.LANCZOS), (x, img_y))
            draw.rectangle([x - 1, img_y - 1, x + w, img_y + h], outline=(200, 200, 200))
            x += w + gap
        y += label_h + row_h + gap
    for p in out_pngs:
        canvas.save(p, dpi=(canvas_w // 10, canvas_w // 10))
        print(f"  {p}  {os.path.getsize(p)/1e6:.1f} MB")
    for p in out_tiffs:
        canvas.save(p, dpi=(canvas_w // 10, canvas_w // 10), compression="tiff_lzw")
        print(f"  {p}  {os.path.getsize(p)/1e6:.1f} MB")


# ---- 300 dpi submission versions (panel cache downscaled internally) ----
SUB = os.path.join(RES, "submission", "Figures")
FM = os.path.join(RES, "figures_main")
os.makedirs(SUB, exist_ok=True); os.makedirs(FM, exist_ok=True)
names300 = {"Figure1": ["Figure1.png", "Fig1_atlas.png"],
            "Figure2": ["Figure2.png", "Fig2_chromatin.png"],
            "Figure3": ["Figure3.png", "Fig3_transcriptome.png"],
            "Figure4": ["Figure4.png", "Fig4_coupling.png"],
            "Figure5": ["Figure5.png", "Fig5_validation_model.png"]}
for base, rows in FIGS.items():
    outs = [os.path.join(SUB, names300[base][0]), os.path.join(FM, names300[base][1])]
    compose(rows, base, 3000, 70, 40, 90, make_font(72), outs, [])

# ---- 600 dpi backup versions ----
OUT6 = os.path.join(RES, "submission", "Figures_600dpi")
os.makedirs(OUT6, exist_ok=True)
for base, rows in FIGS.items():
    compose(rows, base, 6000, 140, 80, 180, make_font(144),
            [os.path.join(OUT6, base + ".png")],
            [os.path.join(OUT6, base + ".tiff")])
print("ALL DONE")
