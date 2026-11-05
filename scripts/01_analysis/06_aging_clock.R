# ============================================================
# Phase 2A - Step 3: chromatin aging clock
# (b) aging TF motifs: donor pseudobulk chromVAR deviations
# (c) per-nucleus chromatin aging score across cell types
# Inputs : chromvar_deviations.rds, cell_metadata.csv
# Outputs: aging_results/ CSVs + plots
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(ggplot2)})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953")
out.dir <- "aging_results"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

dev <- readRDS("chromvar_deviations.rds")
meta <- read.csv("cell_metadata.csv", stringsAsFactors = FALSE)
m <- match(colnames(dev), meta$barcode)
stopifnot(!any(is.na(m)))
meta <- meta[m, ]
cat("deviations:", nrow(dev), "motifs x", ncol(dev), "cells\n"); flush.console()

# ---------- (b) aging TF motifs (donor pseudobulk deviations) ----------
donors <- unique(meta$sample)
group.of <- function(d) unique(meta$group[meta$sample == d])
gd <- factor(sapply(donors, group.of), levels = c("Young", "Aged"))
dvd <- sapply(donors, function(d) Matrix::rowMeans(dev[, meta$sample == d, drop = FALSE]))

motif.res <- do.call(rbind, lapply(seq_len(nrow(dvd)), function(i) {
  x <- dvd[i, ]
  tt <- t.test(x ~ gd)
  data.frame(motif = rownames(dvd)[i],
             meanYoung = mean(x[gd == "Young"]), meanAged = mean(x[gd == "Aged"]),
             diff = diff(tt$estimate), p = tt$p.value)
}))
motif.res$padj <- p.adjust(motif.res$p, method = "BH")
motif.res <- motif.res[order(motif.res$padj), ]
write.csv(motif.res, file.path(out.dir, "motif_activity_Aged_vs_Young_donor_pseudobulk.csv"),
          row.names = FALSE)
aging.motifs <- motif.res$motif[motif.res$padj < 0.05 & motif.res$diff > 0]
youth.motifs  <- motif.res$motif[motif.res$padj < 0.05 & motif.res$diff < 0]
cat("aging motifs (up in Aged):", length(aging.motifs), "| youth motifs:", length(youth.motifs), "\n")
cat("top aging motifs:", paste(head(gsub("^\\S+ ", "", aging.motifs), 12), collapse = ", "), "\n")
flush.console()

# ---------- (c) per-nucleus chromatin aging score ----------
scale.dev <- t(scale(t(dev)))                       # z-score motifs across nuclei
k <- min(50, length(aging.motifs), length(youth.motifs))
sig <- c(head(aging.motifs, k), head(youth.motifs, k))
w <- c(rep(+1, k), rep(-1, k))
score <- as.numeric(colMeans(scale.dev[sig, , drop = FALSE] * w))
meta$agingScore <- score

sc.summ <- aggregate(agingScore ~ celltype + group, data = meta, FUN = mean)
sc.wide <- reshape(sc.summ, idvar = "celltype", timevar = "group", direction = "wide")
colnames(sc.wide) <- c("celltype", "Young", "Aged")
sc.wide$delta <- sc.wide$Aged - sc.wide$Young
sc.wide <- sc.wide[order(-sc.wide$delta), ]
write.csv(sc.wide, file.path(out.dir, "chromatin_aging_score_by_celltype.csv"), row.names = FALSE)
cat("\nchromatin aging score (cell-type ranking, Aged-Young delta):\n")
print(sc.wide); flush.console()

# donor-level paired significance per cell type
ct.test <- do.call(rbind, lapply(sort(unique(meta$celltype)), function(ct) {
  sub <- meta[meta$celltype == ct, ]
  dn <- unique(sub$sample)
  ds <- sapply(dn, function(d) mean(sub$agingScore[sub$sample == d]))
  dg <- factor(sapply(dn, group.of), levels = c("Young", "Aged"))
  tt <- t.test(ds ~ dg)
  data.frame(celltype = ct, delta = diff(tt$estimate), p = tt$p.value)
}))
ct.test$padj <- p.adjust(ct.test$p, "BH")
write.csv(ct.test, file.path(out.dir, "aging_score_donor_test_by_celltype.csv"), row.names = FALSE)
cat("\ndonor-level tests:\n"); print(ct.test); flush.console()

# ---------- plots ----------
meta$celltype <- factor(meta$celltype, levels = sc.wide$celltype)
p <- ggplot(meta[order(meta$celltype), ], aes(celltype, agingScore, fill = group)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = c(Young = "#534AB7", Aged = "#D85A30")) +
  labs(title = "Chromatin aging score by cell type (GSE268953)",
       subtitle = "z-scored chromVAR deviations: aging motifs (+) vs youth motifs (-); donor-level tested",
       x = NULL, y = "chromatin aging score") +
  theme_classic(12) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(out.dir, "aging_score_boxplot.png"), p, width = 10, height = 6, dpi = 150)

if (sum(motif.res$padj < 0.05) >= 5) {
  topm <- motif.res$motif[1:min(15, sum(motif.res$padj < 0.05))]
  vd <- t(scale(t(dvd[match(topm, rownames(dvd)), , drop = FALSE])))
  lab <- gsub("^\\S+ ", "", rownames(vd))
  rownames(vd) <- lab
  png(file.path(out.dir, "top_motifs_donor_heatmap.png"), width = 1800, height = 1400, res = 170)
  pheatmap::pheatmap(vd, annotation_col = data.frame(Group = gd, row.names = colnames(vd)),
                     annotation_colors = list(Group = c(Young = "#534AB7", Aged = "#D85A30")),
                     main = "Top aging TF motifs - donor pseudobulk", fontsize = 9)
  dev.off()
}
cat("PHASE2A_DONE\n")
