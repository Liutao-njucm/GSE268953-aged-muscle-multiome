# ============================================================
# 18b: per-donor FAP expression of clock genes (robustness plot)
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(data.table); library(ggplot2)})
deg.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953/deg"
fig.dir <- "D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953/figures"
d <- fread(file.path(deg.dir, "..", "FAP_clock_per_donor.csv"))
d <- d[!is.na(donor)]
d[, gene := factor(gene, levels = c("NR1D1","NR1D2","PER1","PER3","DBP","TEF","HLF","CRY2",
                                    "ARNTL","NPAS2","IL6","CCL2","LUM","COL1A1","MUSK","PDGFRA"))]
d[, lab := ifelse(grepl("^Young", donor), "Young (20-27 y)", "Aged (60-77 y)")]
p <- ggplot(d, aes(lab, log2(CPM10k + 0.05), colour = lab)) +
  geom_point(position = position_jitter(width = 0.12), size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "black") +
  facet_wrap(~ gene, scales = "free_y", ncol = 8) +
  scale_colour_manual(values = c(`Young (20-27 y)` = "#2E86C1", `Aged (60-77 y)` = "#C0392B"), name = NULL) +
  labs(x = NULL, y = "log2 expression (CPM/10k, FAP pseudobulk)",
       title = "Per-donor expression in FAP: clock repression is donor-consistent") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(fig.dir, "FAP_clock_per_donor.png"), p, width = 10.5, height = 4.4, dpi = 300)
cat("saved FAP_clock_per_donor.png\n")
