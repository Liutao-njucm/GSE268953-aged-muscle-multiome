# ============================================================
# 07: chromatin aging clock analysis
#   a) donor-level motif activity: Aged vs Young (5v5)
#   b) define aging TF signature from global contrast
#   c) per-cell chromatin aging score -> per cell type aging magnitude
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Matrix); library(ggplot2); library(matrixStats)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")
fig.dir <- file.path(out.dir, "figures"); dir.create(fig.dir, showWarnings = FALSE)

dat <- readRDS(file.path(out.dir, "motif_activity_merged.rds"))
z <- dat$z; meta <- dat$meta; motifs <- dat$motif
cat("loaded:", nrow(z), "motifs x", ncol(z), "cells\n")

meta$sample <- factor(meta$sample, levels = c(paste0("Young_", 1:5), paste0("Aged_", 1:5)))
grp_of <- function(s) ifelse(grepl("^Aged", as.character(s)), "Aged", "Young")

# ---- helper: donor-level means of z for a cell subset ----
donor_means <- function(cells) {
  smp <- droplevels(meta$sample[cells])
  D <- sparseMatrix(i = seq_along(cells), j = as.integer(smp), x = 1,
                    dims = c(length(cells), nlevels(smp)))
  cnt <- Matrix::colSums(D)
  M <- as.matrix(z[, cells, drop = FALSE] %*% D)
  M <- sweep(M, 2, cnt, "/")
  colnames(M) <- levels(smp)
  M
}

# ---- helper: 5v5 moderated t-test on a motifs x donors matrix ----
ttest_mat <- function(M) {
  g <- grp_of(colnames(M))
  a <- M[, g == "Aged", drop = FALSE]; y <- M[, g == "Young", drop = FALSE]
  stopifnot(ncol(a) >= 3, ncol(y) >= 3)
  na <- rowSums(!is.na(a)); ny <- rowSums(!is.na(y))
  mA <- rowMeans(a, na.rm = TRUE); mY <- rowMeans(y, na.rm = TRUE)
  vA <- rowVars(a, na.rm = TRUE); vY <- rowVars(y, na.rm = TRUE)
  sp <- sqrt(((na - 1) * vA + (ny - 1) * vY) / (na + ny - 2))
  t <- (mA - mY) / (sp + 1e-12)
  df <- na + ny - 2
  p <- 2 * pt(-abs(t), df)
  p[!is.finite(p)] <- NA
  ok <- !is.na(p) & !is.na(mA) & !is.na(mY)
  out <- data.frame(motif = rownames(M), delta = mA - mY, mean_Aged = mA, mean_Young = mY,
                    t = as.numeric(t), p = as.numeric(p),
                    fdr = NA_real_, stringsAsFactors = FALSE)
  out$fdr[ok] <- p.adjust(out$p[ok], "BH")
  out[ok, , drop = FALSE]
}

# ---- a1) global (all cells) ----
Mg <- donor_means(seq_len(ncol(z)))
res.global <- ttest_mat(Mg)
res.global <- res.global[order(res.global$p), ]
write.csv(res.global, file.path(out.dir, "motif_aging_global_5v5.csv"), row.names = FALSE)
cat("global: FDR<0.05 =", sum(res.global$fdr < 0.05, na.rm = TRUE), "/", nrow(res.global),
    "| tested motifs:", nrow(res.global), "\n")

# ---- a2) per cell type ----
cts <- sort(unique(meta$celltype))
res.ct <- list(); celltype_summary <- list()
for (ct in cts) {
  cells <- which(meta$celltype == ct)
  per_donor <- table(meta$sample[cells])
  if (sum(per_donor[grepl("^Young", names(per_donor))] >= 30) < 3 ||
      sum(per_donor[grepl("^Aged", names(per_donor))] >= 30) < 3) {
    cat("  skip", ct, "(too few donors with >=30 cells)\n"); next
  }
  M <- donor_means(cells)
  r <- ttest_mat(M); r$celltype <- ct
  res.ct[[ct]] <- r
  cat("  ", ct, ": n =", length(cells), "| FDR<0.05 =", sum(r$fdr < 0.05, na.rm = TRUE), "\n")
}
res.ct.all <- do.call(rbind, res.ct)
write.csv(res.ct.all, file.path(out.dir, "motif_aging_by_celltype_5v5.csv"), row.names = FALSE)

# ---- b) aging TF signature ----
sig.file <- file.path(out.dir, "aging_TF_signature.csv")
if (file.exists(sig.file)) {
  sg <- read.csv(sig.file, stringsAsFactors = FALSE)
  up <- intersect(sg$motif[sg$direction == "up_in_Aged"], rownames(z))
  dn <- intersect(sg$motif[sg$direction == "down_in_Aged"], rownames(z))
  cat("signature from motif enrichment of DA peaks:", length(up), "up /", length(dn), "down\n")
} else {
  res.global <- res.global[!is.na(res.global$p), ]
  sig <- res.global[!is.na(res.global$fdr) & res.global$fdr < 0.05, ]
  if (nrow(sig) < 20) {
    cat("  few significant motifs (", nrow(sig), "), fall back to top 100 by p\n")
    sig <- head(res.global[order(res.global$p), ], 100)
  }
  su <- sig[sig$delta > 0, ]; su <- su[order(-su$delta), ]
  sd_ <- sig[sig$delta < 0, ]; sd_ <- sd_[order(sd_$delta), ]
  up <- intersect(head(su$motif, 40), rownames(z))
  dn <- intersect(head(sd_$motif, 40), rownames(z))
  cat("aging TF signature (from motif activity):up =", length(up), "down =", length(dn), "\n")
}
write.csv(data.frame(motif = c(up, dn), direction = c(rep("up_in_Aged", length(up)),
                                                       rep("down_in_Aged", length(dn)))),
          file.path(out.dir, "aging_TF_signature.csv"), row.names = FALSE)

# ---- c) per-cell chromatin aging score ----
zup <- z[up, , drop = FALSE]; zdn <- z[dn, , drop = FALSE]
score <- colMeans(zup) - colMeans(zdn)
score <- (score - mean(score)) / sd(score)
meta$aging_score <- score
write.csv(meta[, c("barcode", "sample", "group", "celltype", "aging_score")],
          file.path(out.dir, "cell_aging_score.csv"), row.names = FALSE)

# per celltype x donor mean score -> 5v5 test
cs <- meta[, c("celltype", "sample", "aging_score")]
agg <- aggregate(aging_score ~ celltype + sample, data = cs, FUN = mean)
n_by <- aggregate(list(n = rep(1L, nrow(cs))),
                  by = list(celltype = cs$celltype, sample = cs$sample), FUN = sum)
agg <- merge(agg, n_by, by = c("celltype", "sample"))
agg <- agg[agg$n >= 30, ]
agg$group <- ifelse(grepl("^Aged", agg$sample), "Aged", "Young")

st <- do.call(rbind, lapply(split(agg, agg$celltype), function(d) {
  if (sum(d$group == "Aged") < 3 || sum(d$group == "Young") < 3) return(NULL)
  tt <- t.test(aging_score ~ group, data = d)
  data.frame(celltype = d$celltype[1],
             mean_Young = mean(d$aging_score[d$group == "Young"]),
             mean_Aged = mean(d$aging_score[d$group == "Aged"]),
             delta = mean(d$aging_score[d$group == "Aged"]) - mean(d$aging_score[d$group == "Young"]),
             p = tt$p.value, n_cells = sum(d$n))
}))
st$fdr <- p.adjust(st$p, "BH")
st <- st[order(-st$delta), ]
write.csv(st, file.path(out.dir, "chromatin_aging_score_by_celltype.csv"), row.names = FALSE)
cat("\n=== chromatin aging score by cell type ===\n"); print(st)

# ---- figures ----
p1 <- ggplot(agg, aes(x = reorder(celltype, aging_score, FUN = mean),
                      y = aging_score, fill = group)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c(Young = "#4C9BE8", Aged = "#E8604C")) +
  labs(x = NULL, y = "chromatin aging score (donor mean)", title = "ATAC-based chromatin aging score") +
  theme_bw(base_size = 12) + coord_flip()
ggsave(file.path(fig.dir, "chromatin_aging_score_by_celltype.png"), p1, width = 8, height = 5.5, dpi = 150)

p2 <- ggplot(meta, aes(x = reorder(celltype, aging_score, FUN = median),
                       y = aging_score, fill = group)) +
  geom_boxplot(outlier.size = 0.2, alpha = 0.8) +
  scale_fill_manual(values = c(Young = "#4C9BE8", Aged = "#E8604C")) +
  labs(x = NULL, y = "chromatin aging score (per cell)", title = "Single-nucleus chromatin aging score") +
  theme_bw(base_size = 12) + coord_flip()
ggsave(file.path(fig.dir, "aging_score_percell_by_celltype.png"), p2, width = 8, height = 5.5, dpi = 150)

# heatmap: aging signature TFs x cell type (donor-averaged z)
top <- c(up, dn)
top <- intersect(top, rownames(z))
hm <- do.call(rbind, lapply(cts, function(ct) {
  cells <- which(meta$celltype == ct)
  if (length(cells) < 50) return(NULL)
  data.frame(celltype = ct, motif = top, z = rowMeans(z[top, cells, drop = FALSE]))
}))
p3 <- ggplot(hm, aes(x = celltype, y = factor(motif, levels = rev(top)), fill = z)) +
  geom_tile() + scale_fill_gradient2(low = "#2C6FA6", mid = "white", high = "#B22222") +
  labs(x = NULL, y = NULL, title = "Aging signature TF motif activity across cell types") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig.dir, "aging_TF_heatmap.png"), p3, width = 8.5, height = 8, dpi = 150)

cat("\nfigures written to", fig.dir, "\n")
