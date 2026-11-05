## 19_main_figures.R  —  Phase 3: auxiliary panels + assemble 5 main figures
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({
  library(data.table); library(ggplot2); library(png); library(grid); library(gridExtra)
})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953")
dir.create("figures_main", showWarnings = FALSE)
COL_UP <- "#C0392B"; COL_DOWN <- "#2E86C1"

## ---------- A. DA peak counts per celltype ----------
da <- fread("DA_peaks_Aged_vs_Young_by_celltype.csv")
da.sig <- da[!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1]
da.sig[, dir := ifelse(log2FoldChange > 0, "Up in Aged", "Down in Aged")]
ct.order <- c("MuSC","Myonuclei","SMC/Pericyte","EC","Macrophage","T/NK cells","FAP")
da.sum <- da.sig[, .N, by = .(celltype, dir)]
pDA <- ggplot(da.sum, aes(factor(celltype, levels = ct.order), N, fill = dir)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("Up in Aged" = COL_UP, "Down in Aged" = COL_DOWN), name = NULL) +
  labs(x = NULL, y = "Differentially accessible peaks\n(padj<0.05, |LFC|>1)",
       title = "Chromatin accessibility changes with age") +
  theme_classic(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                            legend.position = "top", plot.title = element_text(face = "bold", size = 11))
ggsave("figures_main/p_DA_counts.png", pDA, width = 6.2, height = 3.4, dpi = 220)

## ---------- B. DEG counts per celltype ----------
ds <- fread("deg/DEG_summary_by_celltype.csv")
ds.long <- melt(ds[, .(celltype, `n_up_Aged`, `n_down_Aged`)],
                id.vars = "celltype", variable.name = "dir", value.name = "N")
ds.long[, dir := ifelse(dir == "n_up_Aged", "Up in Aged", "Down in Aged")]
ds <- ds[, prop_missing := 1 - n_sig / n_genes_tested]
pDEG <- ggplot(ds.long, aes(factor(celltype, levels = ct.order[ct.order %in% celltype]), N, fill = dir)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5), size = 2.9, colour = "white") +
  scale_fill_manual(values = c("Up in Aged" = COL_UP, "Down in Aged" = COL_DOWN), name = NULL) +
  labs(x = NULL, y = "Differentially expressed genes (padj<0.05)",
       title = "Transcriptional changes with age (pseudobulk, 5v5)") +
  theme_classic(11) + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                            legend.position = "top", plot.title = element_text(face = "bold", size = 11))
ggsave("figures_main/p_DEG_counts.png", pDEG, width = 6.2, height = 3.4, dpi = 220)

## ---------- C. Regulatory chain (FAP, significant TFs) ----------
rc <- fread("regulatory_chain_FAP.csv")
rc.sig <- rc[p_binom < 0.05 & DA_direction == "up"][order(median_rna_lfc)]
rc.sig[, TF := factor(TF, levels = TF)]
pCH <- ggplot(rc.sig, aes(TF, frac_rna_up)) +
  geom_col(aes(fill = frac_rna_up), width = 0.72) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  scale_fill_gradient2(low = "#2E86C1", mid = "#F4D03F", high = "#C0392B",
                       midpoint = 0.65, name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Fraction of linked target genes up-regulated in Aged FAP",
       title = "TF -> open enhancer -> target-gene chain (FAP)") +
  theme_classic(11) + theme(legend.position = "none",
                            plot.title = element_text(face = "bold", size = 11))
ggsave("figures_main/p_regchain.png", pCH, width = 4.6, height = 3.4, dpi = 220)

## ---------- D. Model schematic ----------
sch <- data.frame(
  x = c(0.02, 0.36, 0.36, 0.70, 0.36, 0.70, 0.02),
  y = c(0.60, 0.86, 0.34, 0.86, 0.60, 0.34, 0.08),
  w = c(0.30, 0.28, 0.28, 0.28, 0.28, 0.28, 0.96),
  h = c(0.16, 0.14, 0.14, 0.14, 0.14, 0.14, 0.13),
  lab = c("Ageing + disuse\n(motor-neuron decline)",
          "Clock repressor-arm\ncollapse  NR1D1/PER/DBP\n(FAP, padj 3e-9)",
          "Loss of clock licence\n-> AP-1 / CEBP / NFIA\ntake over enhancers",
          "SASP / fibrosis enhancers\nopen  (1,122 DA peaks)",
          "Chromatin primed\n>90% not transcribed\n(OR 1.85-2.55)",
          "FAP senescence &\nfibrotic niche\n(IL6, LUM, POSTN...)",
          "Myonuclear niche: glycolytic / sarcomeric decline ; MUSK up-regulation"),
  fill = c("#E8E8E8", "#D6EAF8", "#FCF3CF", "#FADBD8", "#FADBD8", "#F5B7B1", "#E8F8F5"),
  stringsAsFactors = FALSE)
pMO <- ggplot(sch, aes(x, y)) +
  geom_rect(aes(xmin = x, xmax = x + w, ymin = y - h/2, ymax = y + h/2, fill = fill),
            colour = "grey30", linewidth = 0.35) +
  geom_text(aes(label = lab), size = 2.55, lineheight = 0.85) +
  annotate("segment", x = 0.32, xend = 0.36, y = 0.66, yend = 0.84, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.79, yend = 0.70, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.64, xend = 0.70, y = 0.84, yend = 0.66, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.64, xend = 0.70, y = 0.36, yend = 0.54, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.27, yend = 0.36, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.52, yend = 0.60, linewidth = 0.5, colour = NA) +
  scale_fill_identity() +
  xlim(0, 1.05) + ylim(0, 1) +
  labs(title = "Working model: clock-licensed senescence enhancers drive the aged FAP niche") +
  theme_void(11) + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
ggsave("figures_main/p_model.png", pMO, width = 7.4, height = 3.2, dpi = 220)
cat("aux panels done\n")

## ================= assemble main figures =================
lab <- function(v) grid.text(v, x = 0.012, y = 0.985, just = c("left","top"),
                             gp = gpar(fontsize = 17, fontface = "bold"))
compose <- function(cells, files, out, widths, heights, labels) {
  png(out, width = 3000, height = round(3000 * sum(heights) / sum(widths)), res = 210)
  lay <- grid.layout(nrow = length(heights), ncol = length(widths),
                     widths = unit(widths, "null"), heights = unit(heights, "null"))
  pushViewport(viewport(layout = lay))
  for (i in seq_along(files)) {
    if (is.na(files[i]) || !file.exists(files[i])) {
      cat("WARN missing:", files[i], "\n"); next
    }
    r <- cells[[i]][1]; cl <- cells[[i]][2]
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = cl))
    img <- readPNG(files[i])
    imgA <- nrow(img) / ncol(img)
    cellA <- (heights[r] / sum(heights)) / (widths[cl] / sum(widths))
    cellA <- cellA * (3000 / (3000 * sum(heights) / sum(widths)))  # = heights[r]/widths[cl]
    cellA <- heights[r] / widths[cl]
    if (imgA > cellA) {  # image taller -> fit height
      grid.raster(img, x = 0.5, y = 0.5, height = 1, width = unit(cellA / imgA, "npc"))
    } else {             # image wider -> fit width
      grid.raster(img, x = 0.5, y = 0.5, width = 1, height = unit(imgA / cellA, "npc"))
    }
    grid.text(labels[i], x = 0.015, y = 0.97, just = c("left","top"),
              gp = gpar(fontsize = 20, fontface = "bold"))
    popViewport()
  }
  dev.off(); cat("wrote", out, "\n")
}

F <- "figures"; FM <- "figures_main"
## Fig 1: atlas
compose(cells = list(c(1,1), c(1,2), c(2,1)),
        files = c("umap_celltype.png", "umap_group.png", "figures/celltype_proportion_stats.png"),
        out = "figures_main/Fig1_atlas.png",
        widths = c(1, 1, 2), heights = c(1, 1.05), labels = c("A","B","C"))
## Fig 2: chromatin
compose(cells = list(c(1,1), c(1,2), c(2,1)),
        files = c("figures_main/p_DA_counts.png", "figures/aging_TF_heatmap.png",
                  "figures/chromatin_aging_score_by_celltype.png"),
        out = "figures_main/Fig2_chromatin.png",
        widths = c(1, 1, 2), heights = c(1, 1.05), labels = c("A","B","C"))
## Fig 3: transcriptome
compose(cells = list(c(1,1), c(2,1), c(2,2)),
        files = c("figures_main/p_DEG_counts.png", "figures/pathway_heatmap.png",
                  "figures/FAP_clock_per_donor.png"),
        out = "figures_main/Fig3_transcriptome.png",
        widths = c(1.15, 1), heights = c(1, 1.25), labels = c("A","B","C"))
## Fig 4: coupling & regulatory chain
compose(cells = list(c(1,1), c(1,2), c(2,1), c(2,2)),
        files = c("figures/coupling_FAP.png", "figures/coupling_class_fraction.png",
                  "figures/clock_TF_decoupling_FAP.png", "figures_main/p_regchain.png"),
        out = "figures_main/Fig4_coupling.png",
        widths = c(1.25, 1), heights = c(1.1, 1), labels = c("A","B","C","D"))
## Fig 5: external validation + model
compose(cells = list(c(1,1), c(2,1), c(3,1)),
        files = c("figures/gtex_concordance_by_module.png", "figures/clock_loci_tracks.png",
                  "figures_main/p_model.png"),
        out = "figures_main/Fig5_validation_model.png",
        widths = c(1), heights = c(1, 1, 0.9), labels = c("A","B","C"))
cat("ALL DONE\n")
