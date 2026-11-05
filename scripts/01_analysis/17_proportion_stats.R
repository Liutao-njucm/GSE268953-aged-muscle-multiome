# ============================================================
# 17: cell-type proportion test (5v5 donor-level, arcsine-sqrt)
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2); library(ggrepel)})

out.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953"
fig.dir <- file.path(out.dir, "figures")

# recompute donor-level proportions from annotated object metadata is heavy;
# use cell_aging_score per-cell csv? instead read the per-sample counts saved earlier
## rebuild from cell_metadata.csv (per-cell table)
meta <- fread(file.path(out.dir, "cell_metadata.csv"))
meta[, sample := as.character(sample)]
tab <- meta[, .N, by = .(sample, celltype)]
tot <- meta[, .N, by = sample]
setnames(tot, "N", "n_total")
tab <- merge(tab, tot, by = "sample")
tab[, prop := N / n_total]
tab[, group := ifelse(grepl("^Aged", sample), "Aged", "Young")]

res <- list()
for (ct in unique(tab$celltype)) {
  d <- tab[celltype == ct]
  py <- d[group == "Young", prop]; pa <- d[group == "Aged", prop]
  if (length(py) < 4 || length(pa) < 4) next
  y <- asin(sqrt(pa)); x <- asin(sqrt(py))
  wt <- suppressWarnings(wilcox.test(x, y))
  tt <- t.test(x, y)
  res[[ct]] <- data.table(celltype = ct,
                          mean_young = mean(py), mean_aged = mean(pa),
                          fc = (mean(pa) + 1e-4) / (mean(py) + 1e-4),
                          log2FC = log2((mean(pa) + 1e-4) / (mean(py) + 1e-4)),
                          p_wilcox = wt$p.value, p_ttest = tt$p.value)
}
r <- rbindlist(res)[, padj := p.adjust(p_ttest, "BH")][order(padj)]
setorder(r, -log2FC)
fwrite(r, file.path(out.dir, "celltype_proportion_stats.csv"))
cat("=== donor-level proportion test (5v5) ===\n")
print(as.data.frame(r), digits = 3)

## figure
dsum <- tab[, .(mean = mean(prop), sd = sd(prop)), by = .(celltype, group)]
lv <- r[order(-log2FC), celltype]
dsum[, celltype := factor(celltype, levels = lv)]
p <- ggplot(dsum, aes(celltype, mean * 100, fill = group)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = (mean - sd) * 100, ymax = (mean + sd) * 100),
                position = position_dodge(0.75), width = 0.25) +
  scale_fill_manual(values = c(Young = "#2E86C1", Aged = "#C0392B"), name = NULL) +
  labs(x = NULL, y = "proportion (%)",
       title = "Cell-type composition shift (donor mean ± SD)",
       subtitle = "Wilcoxon 5v5: FAP p = 0.008, EC p = 0.008") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
ggsave(file.path(fig.dir, "celltype_proportion_stats.png"), p, width = 7.4, height = 4.6, dpi = 300)
cat("saved celltype_proportion_stats.png\n")
