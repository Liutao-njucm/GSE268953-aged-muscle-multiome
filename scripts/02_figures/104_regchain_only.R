## 104_regchain_only.R — 只重绘 Fig 4D（regulatory chain）面板
## 与 48_fix_panels.R 的第 5 节一致，但副标题改为双检验口径
## （正文引用 Fisher padj；14 个显著 TF 同时通过 binomial 与 Fisher）
## 切勿整跑 48：其 5C 节会用旧版轨迹图覆盖 55_redraw_5c.R 的成果。
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2)})
res     <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
pan.dir <- file.path(res, "panels_600dpi"); dir.create(pan.dir, showWarnings = FALSE)
DPI <- 600

rc <- fread(file.path(res, "regulatory_chain_FAP.csv"))
rc.sig <- rc[p_binom < 0.05 & DA_direction == "up"][order(median_rna_lfc)]
fam <- function(t) {
  ifelse(t %in% c("NR1D1","NR1D2","NR2F1","NR2F2","DBP","TEF","HLF","NFIL3"), "Clock",
  ifelse(t %in% c("FOS","JUN","JUNB","JUND","FOSL1","FOSL2","BATF","MAFK","MAFB","MAF"), "AP-1",
  ifelse(t %in% c("CEBPA","CEBPB","CEBPD","CEBPE"), "CEBP", "Other")))
}
rc.sig[, family := fam(TF)]
rc.sig[, TF := factor(TF, levels = TF)]
rc.sig[, nlab := paste0("n=", n_linked_genes)]
fpal <- c(Clock = "#8E44AD", "AP-1" = "#C0392B", CEBP = "#E67E22", Other = "#5D6D7E")
pCH <- ggplot(rc.sig, aes(TF, frac_rna_up, fill = family)) +
  geom_col(width = 0.72, colour = "grey25", linewidth = 0.15) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_text(aes(label = nlab), hjust = -0.12, size = 2.3, colour = "grey30") +
  scale_fill_manual(values = fpal, name = "TF family") +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Fraction of linked target genes up-regulated in Aged FAP",
       title = "TF -> open enhancer -> target-gene chain (FAP)",
       subtitle = "Significant chains (exact binomial p < 0.05;\nall also pass Fisher's exact test);\nn = linked genes; dashed line = 49.3% background") +
  theme_classic(11) +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 11))
ggsave(file.path(pan.dir, "p_regchain_600.png"), pCH, width = 5.4, height = 4.6, dpi = DPI)
cat("4D redrawn:", nrow(rc.sig), "TFs\n")
