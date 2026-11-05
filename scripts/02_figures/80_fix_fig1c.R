# 80_fix_fig1c.R — 重绘 Fig. 1C 面板：修正副标题中错误的 p 值
# 背景：旧副标题硬编码 "FAP p = 0.008, EC p = 0.008"，与实际供体级检验不符
# 实际（脚本 17 输出 celltype_proportion_stats.csv）：FAP Wilcoxon p = 0.095；EC 0.0079；SMC 0.0079；T/NK 0.0159
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
res  <- file.path(proj, "results/GSE268953")
pan  <- file.path(res, "panels_600dpi")
fig  <- file.path(res, "figures")

meta <- fread(file.path(res, "cell_metadata.csv"))
meta[, sample := as.character(sample)]
tab <- meta[, .N, by = .(sample, celltype)]
tot <- meta[, .N, by = sample]; setnames(tot, "N", "n_total")
tab <- merge(tab, tot, by = "sample")
tab[, prop := N / n_total]
tab[, group := ifelse(grepl("^Aged", sample), "Aged", "Young")]

r <- fread(file.path(res, "celltype_proportion_stats.csv"))
r[, padj_wilcox := p.adjust(p_wilcox, "BH")]
fwrite(r, file.path(res, "celltype_proportion_stats.csv"))
cat("=== donor-level proportion tests ===\n")
print(r[order(-log2FC), .(celltype, mean_young, mean_aged, p_wilcox, padj_wilcox, p_ttest, padj)],
      digits = 3)

sub_txt <- sprintf("Wilcoxon (5 vs 5): EC P = %.3f, SMC/pericyte P = %.3f; FAP P = %.3f (ns)",
                   r[celltype == "EC", p_wilcox],
                   r[celltype == "SMC/Pericyte", p_wilcox],
                   r[celltype == "FAP", p_wilcox])
cat("subtitle:", sub_txt, "\n")

dsum <- tab[, .(mean = mean(prop), sd = sd(prop)), by = .(celltype, group)]
lv <- r[order(-log2FC), celltype]
## 其余未纳入检验的细胞类型（每供体核数过少）排在最后，避免 ggplot 把它们折成 "NA" 条
extra <- setdiff(unique(dsum$celltype), lv)
dsum[, celltype := factor(celltype, levels = c(lv, extra))]

p <- ggplot(dsum, aes(celltype, mean * 100, fill = group)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = (mean - sd) * 100, ymax = (mean + sd) * 100),
                position = position_dodge(0.75), width = 0.25) +
  scale_fill_manual(values = c(Young = "#2E86C1", Aged = "#C0392B"), name = NULL) +
  labs(x = NULL, y = "proportion (%)",
       title = "Cell-type composition shift (donor mean \u00b1 SD)",
       subtitle = sub_txt) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")

ggsave(file.path(fig, "celltype_proportion_stats.png"), p, width = 7.4, height = 4.6, dpi = 300)
ggsave(file.path(pan, "celltype_proportion_stats_600.png"), p, width = 7.4, height = 4.6, dpi = 600)
cat("saved panel at 300 and 600 dpi\n")
