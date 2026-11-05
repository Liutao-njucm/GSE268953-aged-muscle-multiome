# ============================================================
# Phase 2A - Step 2: consensus-peak donor-level DA (Aged vs Young)
# Consensus peaks = reduce() of all per-sample peaks; each sample's
# pseudobulk counts mapped to the consensus peak with max overlap.
# Input : atac_pseudobulk_per_sample.rds
# Output: aging_results/DA_peaks_Aged_vs_Young_donor_pseudobulk.csv
# ============================================================
.libPaths(c("C:/Users/wang/R/win-library/4.5", .libPaths()))
suppressMessages({library(Signac); library(GenomicRanges); library(DESeq2)})
setwd("D:/workbuddy/outputs/2026-10-25-20-23-31/results/GSE268953")
out.dir <- "aging_results"
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

x <- readRDS("atac_pseudobulk_per_sample.rds")
pb.list <- x$pb; gr.list <- x$gr

# ---- consensus peaks ----
all.gr <- do.call(c, unname(gr.list))
cons <- reduce(all.gr, minoverlap = 1L)
cat("total per-sample peaks:", length(all.gr), "| consensus peaks:", length(cons), "\n")
flush.console()

# ---- map each sample's peaks to consensus (max-overlap assignment) ----
donors <- colnames(do.call(cbind, pb.list))
pb.cons <- matrix(0L, nrow = length(cons), ncol = length(donors),
                  dimnames = list(NULL, donors))
for (s in names(pb.list)) {
  ov <- findOverlaps(gr.list[[s]], cons, select = "arbitrary")  # one-to-one
  hit <- subjectHits(ov); q <- queryHits(ov)
  agg <- Matrix::rowsum(pb.list[[s]], group = hit, reorder = FALSE)
  pb.cons[as.integer(rownames(agg)), ] <- pb.cons[as.integer(rownames(agg)), , drop = FALSE] + agg
  cat(s, "mapped:", nrow(agg), "peaks\n"); flush.console()
}

# ---- DESeq2 Aged vs Young (donor-level, n=10) ----
grp <- ifelse(grepl("^Young", donors), "Young", "Aged")
grp <- factor(grp, levels = c("Young", "Aged"))
keep <- Matrix::rowSums(pb.cons >= 10) >= 5
cat("consensus peaks tested:", sum(keep), "\n"); flush.console()
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = pb.cons[keep, ],
  colData = data.frame(group = grp, row.names = donors),
  design = ~ group)
dds <- DESeq2::DESeq(dds, quiet = TRUE)
res <- DESeq2::results(dds, contrast = c("group", "Aged", "Young"))
res.df <- as.data.frame(res)
res.df$peak <- paste0(as.character(seqnames(cons)[keep]), ":",
                      start(cons)[keep], "-", end(cons)[keep])
res.df <- res.df[order(res.df$padj), ]
write.csv(res.df, file.path(out.dir, "DA_peaks_Aged_vs_Young_donor_pseudobulk.csv"),
          row.names = FALSE)
cat("DA peaks FDR<0.05:", sum(res.df$padj < 0.05, na.rm = TRUE),
    "| up in Aged:", sum(res.df$padj < 0.05 & res.df$log2FoldChange > 0, na.rm = TRUE),
    "| down in Aged:", sum(res.df$padj < 0.05 & res.df$log2FoldChange < 0, na.rm = TRUE), "\n")
cat("STEP2_DONE\n")
