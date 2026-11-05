# ============================================================
# 08: differential chromatin accessibility (consensus peaks)
#   donor-level pseudobulk (each sample = 1 donor), DESeq2 5v5
#   run per major cell type to avoid composition confounding
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Matrix); library(DESeq2); library(GenomicRanges)})

proj <- "D:/workbuddy/outputs/2026-10-25-20-23-31"
out.dir <- file.path(proj, "results/GSE268953")
tmp.dir <- file.path(out.dir, "atac_tmp")

samples <- c(paste0("Young_", 1:5), paste0("Aged_", 1:5))
pbs <- lapply(samples, function(s) readRDS(file.path(tmp.dir, paste0("pb_", s, ".rds"))))
names(pbs) <- samples
cons <- readRDS(file.path(out.dir, "consensus_peaks.rds"))
cat("consensus peaks:", length(cons), "\n")

# ---- restrict to regions called in ALL 10 samples ----
# otherwise regions called only in some donors create fake
# "less accessible in Aged" signal
support <- readRDS(file.path(out.dir, "consensus_support.rds"))
sel <- which(support == 10)
cat("regions called in all 10 samples:", length(sel), "/", length(cons), "\n")
cons.all <- cons
cons <- cons[sel]
sel.id <- sel

celltypes <- sort(unique(unlist(lapply(pbs, function(x) colnames(x$pb)))))
grp <- factor(ifelse(grepl("^Aged", samples), "Aged", "Young"), levels = c("Young", "Aged"))

all.res <- list()
for (ct in celltypes) {
  mats <- lapply(pbs, function(x) if (ct %in% colnames(x$pb)) x$pb[, ct] else NULL)
  have <- !sapply(mats, is.null)
  ncell <- sapply(pbs, function(x) if (ct %in% colnames(x$pb)) x$celltype_n[[ct]] else 0)
  ok <- have & ncell >= 30
  if (sum(ok & grp == "Young") < 4 || sum(ok & grp == "Aged") < 4) {
    cat("skip", ct, ": donors with >=30 cells ->", sum(ok), "\n"); next
  }
  M <- do.call(cbind, mats[ok])
  M <- M[sel.id, , drop = FALSE]
  rownames(M) <- as.character(sel.id)
  colnames(M) <- samples[ok]
  g <- droplevels(grp[ok])

  # keep accessible peaks only (reduce DESeq2 burden)
  keep <- Matrix::rowSums(M >= 5) >= 3
  M <- M[keep, , drop = FALSE]
  if (nrow(M) > 120000) {
    top <- order(Matrix::rowSums(M), decreasing = TRUE)[1:120000]
    M <- M[sort(top), , drop = FALSE]
  }
  cat(ct, ": ", nrow(M), "peaks x", ncol(M), "donors\n")

  dds <- DESeqDataSetFromMatrix(round(as.matrix(M)), DataFrame(group = g), ~group)
  dds <- DESeq(dds, quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("group", "Aged", "Young")))
  r$peak <- rownames(r); r$celltype <- ct
  r <- r[!is.na(r$padj), ]
  all.res[[ct]] <- r
  cat("  sig (padj<0.05, |LFC|>1):", sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE),
      "| up in Aged:", sum(r$padj < 0.05 & r$log2FoldChange > 1, na.rm = TRUE), "\n")
  rm(dds, M); invisible(gc())
}

res <- do.call(rbind, all.res)
gr <- as.data.frame(cons.all)[, c("seqnames", "start", "end")]
names(gr) <- c("chr", "start", "end")
gr$peak <- as.character(seq_len(nrow(gr)))
res <- merge(res, gr, by = "peak", all.x = TRUE)
res$peak_id <- paste0(res$chr, "-", res$start, "-", res$end)
res <- res[order(res$celltype, res$padj), ]
write.csv(res, file.path(out.dir, "DA_peaks_Aged_vs_Young_by_celltype.csv"), row.names = FALSE)
cat("saved DA table:", nrow(res), "rows\n")
