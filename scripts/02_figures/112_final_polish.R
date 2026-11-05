## 112_final_polish.R — 全图终审后的三个面板重绘（2026-11-05）
## ① p_DA_counts_600：y 轴标题换两行（原被上下裁掉 "D"）
## ② aging_TF_heatmap_600：~90 个 TF 行标签原 10pt 互相重叠 → 5.2pt
## ③ p_model_600："+ MUSK compensation" → "; MUSK up-regulation"（与正文口径一致，
##    正文 Results/Discussion 均为中性表述，且 compensation 属禁用词族）
## 切勿整跑 48（会覆盖 55 的 v3 轨迹图）。
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2)})
res     <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
pan.dir <- file.path(res, "panels_600dpi")
DPI <- 600
COL_UP <- "#C0392B"; COL_DOWN <- "#2E86C1"

## ---- ① Fig 2A: DA peak counts ----
da <- fread(file.path(res, "DA_peaks_Aged_vs_Young_by_celltype.csv"))
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
ggsave(file.path(pan.dir, "p_DA_counts_600.png"), pDA, width = 6.2, height = 3.4, dpi = DPI)
cat("2A done\n")

## ---- ② Fig 2B: TF motif heatmap（数据来自已保存的 merged rds，不需 chromVAR）----
dat <- readRDS(file.path(res, "motif_activity_merged.rds"))
z <- dat$z; meta <- dat$meta; rm(dat); invisible(gc())
sg <- read.csv(file.path(res, "aging_TF_signature.csv"), stringsAsFactors = FALSE)
up <- intersect(sg$motif[sg$direction == "up_in_Aged"], rownames(z))
dn <- intersect(sg$motif[sg$direction == "down_in_Aged"], rownames(z))
cts <- sort(unique(meta$celltype))
hm <- do.call(rbind, lapply(cts, function(ct) {
  cells <- which(meta$celltype == ct)
  if (length(cells) < 50) return(NULL)
  data.frame(celltype = ct, motif = c(up, dn),
             z = rowMeans(z[c(up, dn), cells, drop = FALSE]))
}))
top <- c(up, dn)
p3 <- ggplot(hm, aes(x = celltype, y = factor(motif, levels = rev(top)), fill = z)) +
  geom_tile() + scale_fill_gradient2(low = "#2C6FA6", mid = "white", high = "#B22222") +
  labs(x = NULL, y = NULL, title = "Aging signature TF motif activity across cell types") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                   axis.text.y = element_text(size = 5.2))
ggsave(file.path(pan.dir, "aging_TF_heatmap_600.png"), p3, width = 8.5, height = 8, dpi = DPI)
rm(z, meta); invisible(gc())
cat("2B done,", nrow(hm), "rows\n")

## ---- ③ Fig 5D: working model（与 48 第 4 节一致，仅底行文字修正）----
sch <- data.frame(
  x = c(0.02, 0.36, 0.36, 0.70, 0.36, 0.70, 0.02),
  y = c(0.60, 0.86, 0.34, 0.86, 0.60, 0.34, 0.08),
  w = c(0.30, 0.28, 0.28, 0.28, 0.28, 0.28, 0.96),
  h = c(0.16, 0.14, 0.14, 0.14, 0.14, 0.14, 0.13),
  lab = c("Ageing + disuse\n(motor-neuron decline)",
          "Clock repressor-arm\ncollapse NR1D1/PER/DBP\n(FAP, padj 3e-9)",
          "Loss of clock licence\n-> AP-1 / CEBP / NFIA\ntake over enhancers",
          "SASP / fibrosis enhancers\nopen (1,122 DA peaks)",
          "Chromatin primed\n>90% not transcribed\n(OR 1.85-2.55)",
          "FAP senescence &\nfibrotic niche\n(IL6, PRG4, ADAMTS16)",
          "Myonuclear niche: glycolytic / sarcomeric decline; MUSK up-regulation"),
  fill = c("#E8E8E8", "#D6EAF8", "#FCF3CF", "#FADBD8", "#FADBD8", "#F5B7B1", "#E8F8F5"),
  stringsAsFactors = FALSE)
pMO <- ggplot(sch, aes(x, y)) +
  geom_rect(aes(xmin = x, xmax = x + w, ymin = y - h/2, ymax = y + h/2, fill = fill),
            colour = "grey30", linewidth = 0.35) +
  geom_text(aes(x = x + w/2, y = y, label = lab), size = 2.45, lineheight = 0.85) +
  annotate("segment", x = 0.32, xend = 0.36, y = 0.66, yend = 0.84, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.79, yend = 0.70, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.64, xend = 0.70, y = 0.84, yend = 0.66, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.64, xend = 0.70, y = 0.36, yend = 0.54, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.53, yend = 0.41, linewidth = 0.5,
           arrow = arrow(length = unit(2.2, "mm"))) +
  scale_fill_identity() + xlim(0, 1.02) + ylim(0, 1) +
  labs(title = "Working model: clock-licensed senescence enhancers drive the aged FAP niche") +
  theme_void(11) + theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
ggsave(file.path(pan.dir, "p_model_600.png"), pMO, width = 7.4, height = 3.2, dpi = DPI)
cat("5D done\n")
